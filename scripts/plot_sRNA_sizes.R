#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(cowplot)

sRNA.Allsense = read_tsv("results/Botrytis_sense_size_profile.tsv") %>% filter(SIZE >= 18)
sRNA.Allantisense = read_tsv("results/Botrytis_antisense_size_profile.tsv") %>% filter(SIZE >= 18)

sRNA.sense <- sRNA.Allsense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "TE.Unknown" & FeatureType != "TE.Simple_repeat" & FeatureType != "TE.RC" & FeatureType != "TE.Low_complexity" &
            FeatureType != "TE.rRNA" & FeatureType != "TE.Satellite" & FeatureType != "TE.Retroposon")
sRNA.noExon <- sRNA.sense %>% filter(FeatureType != 'exon' & FeatureType != "None" & FeatureType != "intron")

negate <- function(x) ( -1.0 * x)
add_antisense <- function(x) ( paste0(x,".antisense"))
sRNA.antisense <- sRNA.Allantisense %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "FeatureType", values_to = "count") %>%
    filter(FeatureType != "TE.Unknown" & FeatureType != "TE.Simple_repeat" & FeatureType != "TE.RC" & FeatureType != "TE.Low_complexity" &
           FeatureType != "TE.rRNA" & FeatureType != "TE.Satellite" & FeatureType != "TE.Retroposon") %>%
    mutate_at(vars(-c("SIZE","FeatureType")),negate) %>% mutate_at(c("FeatureType"),add_antisense)
sRNA.antinoExon <- sRNA.antisense %>% filter(FeatureType != 'exon.antisense' & FeatureType != "None.antisense" & FeatureType != "intron.antisense")

sRNA.combined <- bind_rows(sRNA.antinoExon,sRNA.noExon) %>% arrange(desc(FeatureType))

pairedcolors = brewer.pal(12, "Paired")
green1 = pairedcolors[3]
green2 = pairedcolors[4]
pairedcolors[3] = "#FFCC00"
pairedcolors[4] = "#B15928"
pairedcolors[11] = green1
pairedcolors[12] = green2

pN<-ggplot(sRNA.combined,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA.combined$SIZE)) +
    theme_cowplot(12) + scale_fill_manual(values = pairedcolors)
pN

sRNA.Exon <- bind_rows(sRNA.sense %>% filter(FeatureType == 'exon' | FeatureType == "intron"),
                       sRNA.antisense %>% filter(FeatureType == 'exon.antisense' | FeatureType == "intron.antisense"))

pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_brewer(palette="Paired") +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined_NoExon$SIZE))
pE

sRNA.Other <- sRNA.sense %>% filter(FeatureType == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=FeatureType)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_brewer(palette="Set2") +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA.combined_NoExon$SIZE))
pO


pg <- plot_grid(pN, pE, pO, labels = c('NonCoding', 'Coding', "Intergenic"), label_size = 12,ncol = 1,align = "v")
ggsave("size_profile_plots.pdf",height=18,width=12)

# pie chart time

ggplot(data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)
