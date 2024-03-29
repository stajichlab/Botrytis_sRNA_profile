#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(cowplot)

# some definitions
pairedcolors = brewer.pal(12, "Paired")
green1 = pairedcolors[3]
green2 = pairedcolors[4]
pairedcolors[3] = "#FFCC00"
pairedcolors[4] = "#B15928"
pairedcolors[9] = green1
pairedcolors[10] = green2
pairedcolors[13] = '#FFFF11'
pairedcolors[14] = '#AA0011'

exonColors = brewer.pal(6, "Set1")
exonColors[3] = exonColors[5]
IntergenicColors = rev(brewer.pal(3,"Greys"))

piecolors = pairedcolors
piecolors[1] = exonColors[1]
piecolors[2] = exonColors[2]
piecolors[3] = exonColors[3]
piecolors[4] = IntergenicColors[1]
piecolors[5] = pairedcolors[1]
piecolors[6] = pairedcolors[3]
piecolors[7] = pairedcolors[5]
piecolors[8] = pairedcolors[7]
piecolors[9] = pairedcolors[9]
piecolors[10] = pairedcolors[11]

print(exonColors[3])

negate <- function(x) ( -1.0 * x)
add_antisense <- function(x) ( paste0(x,".antisense"))
# setup default plot out
pdf("plots/Combo_all_plots_uniq.pdf")

# read in the "all data"
sRNA.Allsense = read_tsv("results/Botrytis_sense_size_profile_uniq.tsv") %>% filter(SIZE >= 18)
sRNA.Allantisense = read_tsv("results/Botrytis_antisense_size_profile_uniq.tsv") %>% filter(SIZE >= 18)

sRNA.sense <- sRNA.Allsense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count")
sRNA.noExon <- sRNA.sense %>% filter(FeatureType != 'exon' & FeatureType != "None" & FeatureType != "intron")


sRNA.antisense <- sRNA.Allantisense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    mutate_at(vars(-c("SIZE","FeatureType")),negate) %>% mutate_at(c("FeatureType"),add_antisense)
sRNA.antinoExon <- sRNA.antisense %>% filter(FeatureType != 'exon.antisense' & FeatureType != "None.antisense" & FeatureType != "intron.antisense")

sRNA.combined <- bind_rows(sRNA.antinoExon,sRNA.noExon) %>% arrange(desc(FeatureType))



pN<-ggplot(sRNA.combined,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    theme_cowplot(12) + scale_fill_manual(values = pairedcolors) + ggtitle("Combined all Libraries sRNA profile NonCoding Features for Unique Mapping")
pN

sRNA.Exon <- bind_rows(sRNA.sense %>% filter(FeatureType == 'exon' | FeatureType == "intron"),
                       sRNA.antisense %>% filter(FeatureType == 'exon.antisense' | FeatureType == "intron.antisense"))

pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = exonColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Combined all Libraries sRNA profile Exon,Intron Features for Unique Mapping")
pE

sRNA.Other <- sRNA.sense %>% filter(FeatureType == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = IntergenicColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Combined all Libraries sRNA profile for Non-Genic location for Unique Mapping")
pO


pg <- plot_grid(pN, pE, pO, labels = c('A', 'B', "C"), label_size = 12,ncol = 1,align = "v")
ggsave("plots/All_Libraries_size_profile_uniq.pdf",pg,height=18,width=12)

# pie chart time
# make antisense counts, collapse sRNA read length

sRNA.antisenseWoExon <- sRNA.Allsense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "exon") %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count))

exonAntisense <- sRNA.antisense %>% filter(FeatureType == "exon.antisense") %>% group_by(FeatureType) %>%
    mutate_at(vars(c(count)),negate) %>% summarize(countAllSize = sum(count))
# combine all sense feature with antisense features, but leave exons.antisense out as separate category
sRNA.combinedAll <- bind_rows(sRNA.sense %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count)),
                              exonAntisense,sRNA.antisenseWoExon) %>% group_by(FeatureType) %>% summarize(countAll = sum(countAllSize))

Total = sum(sRNA.combinedAll %>% select(countAll))

sRNA.combinedAll <- sRNA.combinedAll %>% mutate(Perc=100*(countAll/Total)) # add percentage in
pie <- ggplot(sRNA.combinedAll, aes(x="", y=countAll, fill=FeatureType)) +
    geom_bar(stat="identity", width=1) + geom_text(aes(label = paste0(sprintf("%.1f",Perc), "%")), position = position_stack(vjust=0.5), size=3) +
    coord_polar("y", start=0) +  scale_fill_manual(values = piecolors) + ggtitle("Total for all Libraries sRNA profile for Unique Mapping") +
    theme_classic() +  labs(x = NULL, y = NULL, fill = NULL) +   theme(axis.line = element_blank(),
                                                                       axis.text = element_blank(),
                                                                       axis.ticks = element_blank())

pie
ggsave("plots/All_Libraries_Piechart_uniq.pdf",pie)

######

# B
sRNA.Bsense = read_tsv("results/Botrytis_sense_ExprB_size_profile_uniq.tsv") %>% filter(SIZE >= 18)
sRNA.Bantisense = read_tsv("results/Botrytis_antisense_ExprB_size_profile_uniq.tsv") %>% filter(SIZE >= 18)

sRNA.sense <- sRNA.Bsense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count")
sRNA.noExon <- sRNA.sense %>% filter(FeatureType != 'exon' & FeatureType != "None" & FeatureType != "intron")

sRNA.antisense <- sRNA.Bantisense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    mutate_at(vars(-c("SIZE","FeatureType")),negate) %>% mutate_at(c("FeatureType"),add_antisense)
sRNA.antinoExon <- sRNA.antisense %>% filter(FeatureType != 'exon.antisense' & FeatureType != "None.antisense" & FeatureType != "intron.antisense")

sRNA.combined <- bind_rows(sRNA.antinoExon,sRNA.noExon) %>% arrange(desc(FeatureType))

pN<-ggplot(sRNA.combined,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    theme_cowplot(12) + scale_fill_manual(values = pairedcolors) + ggtitle("Botrytis Mycelium sRNA profile NonCoding Feature for Unique Mapping")
pN

sRNA.Exon <- bind_rows(sRNA.sense %>% filter(FeatureType == 'exon' | FeatureType == "intron"),
                       sRNA.antisense %>% filter(FeatureType == 'exon.antisense' | FeatureType == "intron.antisense"))

pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = exonColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Botrytis Mycelium Libraries sRNA profile Exon,Intron Features for Unique Mapping")
pE

sRNA.Other <- sRNA.sense %>% filter(FeatureType == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = IntergenicColors)  +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Botrytis Mycelium Libraries sRNA profile for Nongenic location for Unique Mapping")
pO

pg <- plot_grid(pN, pE, pO, labels = c('A', 'B', "C"), label_size = 12,ncol = 1,align = "v")
ggsave("plots/Mycelium_size_profile_uniq.pdf",pg,height=18,width=12)

# pie chart time
# make antisense counts, collapse sRNA read length

sRNA.antisenseWoExon <- sRNA.Bsense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "exon") %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count))

exonAntisense <- sRNA.antisense %>% filter(FeatureType == "exon.antisense") %>% group_by(FeatureType) %>%
    mutate_at(vars(c(count)),negate) %>% summarize(countAllSize = sum(count))
# combine all sense feature with antisense features, but leave exons.antisense out as separate category
sRNA.combinedAll <- bind_rows(sRNA.sense %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count)),
                              exonAntisense,sRNA.antisenseWoExon) %>% group_by(FeatureType) %>% summarize(countAll = sum(countAllSize))

Total = sum(sRNA.combinedAll %>% select(countAll))

sRNA.combinedAll <- sRNA.combinedAll %>% mutate(Perc=100*(countAll/Total)) # add percentage in
pie <- ggplot(sRNA.combinedAll, aes(x="", y=countAll, fill=FeatureType)) +
    geom_bar(stat="identity", width=1) + geom_text(aes(label = paste0(sprintf("%.1f",Perc), "%")), position = position_stack(vjust=0.5), size=3) +
    coord_polar("y", start=0) +  scale_fill_manual(values = piecolors) + ggtitle("Botrytis Mycelium sRNA profile for Unique Mapping") +
    theme_classic() +  labs(x = NULL, y = NULL, fill = NULL) +   theme(axis.line = element_blank(),
                                                                       axis.text = element_blank(),
                                                                       axis.ticks = element_blank())

pie
ggsave("plots/Mycelium_Piechart_uniq.pdf",pie,height=8,width=8)


## Infection
sRNA.Isense = read_tsv("results/Botrytis_sense_ExprI_size_profile_uniq.tsv") %>% filter(SIZE >= 18)
sRNA.Iantisense = read_tsv("results/Botrytis_antisense_ExprI_size_profile_uniq.tsv") %>% filter(SIZE >= 18)

sRNA.sense <- sRNA.Isense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count")
sRNA.noExon <- sRNA.sense %>% filter(FeatureType != 'exon' & FeatureType != "None" & FeatureType != "intron")

sRNA.antisense <- sRNA.Iantisense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    mutate_at(vars(-c("SIZE","FeatureType")),negate) %>% mutate_at(c("FeatureType"),add_antisense)
sRNA.antinoExon <- sRNA.antisense %>% filter(FeatureType != 'exon.antisense' & FeatureType != "None.antisense" & FeatureType != "intron.antisense")

sRNA.combined <- bind_rows(sRNA.antinoExon,sRNA.noExon) %>% arrange(desc(FeatureType))

pN<-ggplot(sRNA.combined,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    theme_cowplot(12) + scale_fill_manual(values = pairedcolors) + ggtitle("Tomato Infection sRNA profile NonCoding Feature for Unique Mapping")
pN

sRNA.Exon <- bind_rows(sRNA.sense %>% filter(FeatureType == 'exon' | FeatureType == "intron"),
                       sRNA.antisense %>% filter(FeatureType == 'exon.antisense' | FeatureType == "intron.antisense"))

pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = exonColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Tomato Infection Libraries sRNA profile Exon,Intron Features for Unique Mapping")
pE

sRNA.Other <- sRNA.sense %>% filter(FeatureType == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = IntergenicColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Tomato Infection Libraries sRNA profile for Nongenic location for Unique Mapping")
pO

pg <- plot_grid(pN, pE, pO, labels = c('A', 'B', "C"), label_size = 12,ncol = 1,align = "v")
ggsave("plots/TomatoInfection_size_profile_uniq.pdf",pg,height=18,width=12)

# pie chart time
# make antisense counts, collapse sRNA read length

sRNA.antisenseWoExon <- sRNA.Isense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "exon") %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count))

exonAntisense <- sRNA.antisense %>% filter(FeatureType == "exon.antisense") %>% group_by(FeatureType) %>%
    mutate_at(vars(c(count)),negate) %>% summarize(countAllSize = sum(count))
# combine all sense feature with antisense features, but leave exons.antisense out as separate category
sRNA.combinedAll <- bind_rows(sRNA.sense %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count)),
                              exonAntisense,sRNA.antisenseWoExon) %>% group_by(FeatureType) %>% summarize(countAll = sum(countAllSize))

Total = sum(sRNA.combinedAll %>% select(countAll))

sRNA.combinedAll <- sRNA.combinedAll %>% mutate(Perc=100*(countAll/Total)) # add percentage in
pie <- ggplot(sRNA.combinedAll, aes(x="", y=countAll, fill=FeatureType)) +
    geom_bar(stat="identity", width=1) + geom_text(aes(label = paste0(sprintf("%.1f",Perc), "%")), position = position_stack(vjust=0.5), size=3) +
    coord_polar("y", start=0) +  scale_fill_manual(values = piecolors) + ggtitle("Botrytis Tomato Infection sRNA profile for Unique Mapping") +
    theme_classic() +  labs(x = NULL, y = NULL, fill = NULL) +   theme(axis.line = element_blank(),
                                                                       axis.text = element_blank(),
                                                                       axis.ticks = element_blank())

pie
ggsave("plots/TomatoInfection_Piechart_uniq.pdf",pie,width=8,height=8)


####
# Mock
sRNA.Msense = read_tsv("results/Botrytis_sense_ExprM_size_profile_uniq.tsv") %>% filter(SIZE >= 18)
sRNA.Mantisense = read_tsv("results/Botrytis_antisense_ExprM_size_profile_uniq.tsv") %>% filter(SIZE >= 18)

sRNA.sense <- sRNA.Msense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count")
sRNA.noExon <- sRNA.sense %>% filter(FeatureType != 'exon' & FeatureType != "None" & FeatureType != "intron")

sRNA.antisense <- sRNA.Mantisense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    mutate_at(vars(-c("SIZE","FeatureType")),negate) %>% mutate_at(c("FeatureType"),add_antisense)
sRNA.antinoExon <- sRNA.antisense %>% filter(FeatureType != 'exon.antisense' & FeatureType != "None.antisense" & FeatureType != "intron.antisense")

sRNA.combined <- bind_rows(sRNA.antinoExon,sRNA.noExon) %>% arrange(desc(FeatureType))

pN<-ggplot(sRNA.combined,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    theme_cowplot(12) + scale_fill_manual(values = pairedcolors) + ggtitle("Tomato Mock Infection sRNA profile NonCoding Feature for Unique Mapping")
pN

pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = exonColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Tomato Mock Infection Libraries sRNA profile Exon,Intron Features for Unique Mapping")
pE

sRNA.Other <- sRNA.sense %>% filter(FeatureType == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_manual(values = IntergenicColors) +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    ggtitle("Tomato Mock Infection Libraries sRNA profile for Nongenic location for Unique Mapping")
pO

pg <- plot_grid(pN, pE, pO, labels = c('A', 'B', "C"), label_size = 12,ncol = 1,align = "v")
ggsave("plots/MockInfection_size_profile_uniq.pdf",pg,height=18,width=12)

# pie chart time
# make antisense counts, collapse sRNA read length

sRNA.antisenseWoExon <- sRNA.Isense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "exon") %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count))

exonAntisense <- sRNA.antisense %>% filter(FeatureType == "exon.antisense") %>% group_by(FeatureType) %>%
    mutate_at(vars(c(count)),negate) %>% summarize(countAllSize = sum(count))
# combine all sense feature with antisense features, but leave exons.antisense out as separate category
sRNA.combinedAll <- bind_rows(sRNA.sense %>% group_by(FeatureType) %>% summarize(countAllSize = sum(count)),
                              exonAntisense,sRNA.antisenseWoExon) %>% group_by(FeatureType) %>% summarize(countAll = sum(countAllSize))

Total = sum(sRNA.combinedAll %>% select(countAll))

sRNA.combinedAll <- sRNA.combinedAll %>% mutate(Perc=100*(countAll/Total)) # add percentage in
pie <- ggplot(sRNA.combinedAll, aes(x="", y=countAll, fill=FeatureType)) +
    geom_bar(stat="identity", width=1) + geom_text(aes(label = paste0(sprintf("%.1f",Perc), "%")), position = position_stack(vjust=0.5), size=3) +
    coord_polar("y", start=0) +  scale_fill_manual(values = piecolors) + ggtitle("Mock Tomato Infection sRNA profile for Unique Mapping") +
    theme_classic() +  labs(x = NULL, y = NULL, fill = NULL) +   theme(axis.line = element_blank(),
                                                                       axis.text = element_blank(),
                                                                       axis.ticks = element_blank())

pie
ggsave("plots/MockInfection_Piechart_uniq.pdf",pie,width=8,height=8)
