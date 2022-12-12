#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(cowplot)

sRNA.All = read_tsv("results/Botrytis_size_profile.tsv") %>% filter(SIZE >= 18)
sRNA <- sRNA.All %>% select(-c('TOTAL')) %>% pivot_longer(!SIZE, names_to = "Type", values_to = "count")
sRNA.noExon <- sRNA %>% filter(Type != 'exon')
pN<-ggplot(sRNA.noExon,aes(SIZE,count,fill=Type)) +
    geom_bar(stat='identity',position="stack") + scale_x_discrete(limit = unique(sRNA$SIZE)) +
    scale_fill_brewer(palette="Set3") +  theme_cowplot(12)
pN

sRNA.Exon <- sRNA %>% filter(Type == 'exon')
pE<-ggplot(sRNA.Exon,aes(SIZE,count,fill=Type)) +
    geom_bar(stat='identity',position="stack") +
    scale_fill_brewer(palette="Set3") +  theme_cowplot(12)  + scale_x_discrete(limit = unique(sRNA$SIZE)) 
pE

plot_grid(pN, pE, labels = c('A', 'B'), label_size = 12,ncol = 1,align = "v")
