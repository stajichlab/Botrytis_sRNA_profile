#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(cowplot)

sRNA.All = read_tsv("results/Botrytis_size_profile.tsv") %>% filter(SIZE >= 18)
sRNA <- sRNA.All %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "Type", values_to = "count") %>%
    filter(Type != "TE.Unknown" & Type != "TE.Simple_repeat" & Type != "TE.RC" & Type != "TE.Low_complexity")
sRNA.noExon <- sRNA %>% filter(Type != 'exon' & Type != "None")
pN<-ggplot(sRNA.noExon,aes(SIZE,count,fill=Type)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA$SIZE)) +
    scale_fill_brewer(palette="RdYlGn") +  theme_cowplot(12)
pN

sRNA.Exon <- sRNA %>% filter(Type == 'exon')
pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=Type)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_brewer(palette="Paired") +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA$SIZE))
pE

sRNA.Other <- sRNA %>% filter(Type == 'None')
pO<-ggplot(sRNA.Other,aes(SIZE,count,fill=Type)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_brewer(palette="Set2") +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA$SIZE))
pO


pg <- plot_grid(pN, pE, pO, labels = c('NonCoding', 'Coding', "Intergenic"), label_size = 12,ncol = 1,align = "v")
ggsave("size_profile_plots.pdf",height=16,width=12)
