#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(vroom)

sRNA.All = vroom("results/Botrytis.tomato_subtract.tsv.gz") %>% filter(LENGTH >= 18 & LENGTH <= 30) %>%
    select(c(LENGTH,TOTAL_COUNT))
sRNA.count <- sRNA.All %>% group_by(LENGTH) %>% summarise(Count = sum(TOTAL_COUNT))
sRNA.count

p<-ggplot(sRNA.count, aes(x=LENGTH, y=Count)) + geom_bar(stat='identity', width=1,color = "#000000", fill = "#0099F8") +
    theme_cowplot(12)
ggsave("Botrytis_tomatosubtract.pdf",p,width=12)
