#!/usr/bin/env Rscript
library(tidyverse)
library(vroom)

sRNA_Bc <- vroom("results/Botrytis_matchAligned.out.reads.tsv.gz")  %>% filter(LENGTH > 18)
write_tsv(sRNA_Bc %>% filter(TOTAL_COUNT > 100),"results/Botrytis_matchAligned.min500.tsv.gz")

sRNA_To <- vroom("results/Tomato_matchAligned.out.reads.tsv.gz") %>% filter(LENGTH > 18)
write_tsv(sRNA_To %>% filter(TOTAL_COUNT > 500),"results/Tomato_matchAligned.min500.tsv.gz")

sRNA_To_mockAdj <- sRNA_To %>% mutate(SumMockTo = M12As+M16As+M24As,
                                      SumBcOnly = B16As+B16Bs,
                                      SumInfx  = I12As+I12Bs+I12Cs+I16As+I16Bs+I16Cs+I24Bs+I24Ds)

TomNoBcInBcOnly <- sRNA_To_mockAdj %>%
                   filter( SumBcOnly <= 10 & SumMockTo > 0) %>%
                   arrange(desc(TOTAL_COUNT))

write_tsv(TomNoBcInBcOnly %>% filter(TOTAL_COUNT > 100),"results/Tomato_matchAligned.mockBcFilter.min100.tsv.gz")

sRNA_Bc_mockAdj <- sRNA_Bc %>% mutate(SumMockTo = M12As+M16As+M24As,
                                      SumBcOnly = B16As+B16Bs,
                                      SumInfx  = I12As+I12Bs+I12Cs+I16As+I16Bs+I16Cs+I24Bs+I24Ds)

BcNoBcInMock <- sRNA_Bc_mockAdj %>% filter( SumMockTo <= 10 ) %>% arrange(desc(TOTAL_COUNT))

write_tsv(sRNA_Bc_mockAdj %>% filter(SumInfx > 100),"results/Botrytis_matchAligned.Infx100.tsv.gz")
write_tsv(BcNoBcInMock,"results/Botrytis_matchAligned.mockToFilter.tsv.gz")

sRNA_Bc_mockAdj_strict <- sRNA_Bc_mockAdj %>% filter(SumBcOnly > 5 & SumMockTo < 5) %>% arrange(desc(SumInfx))
sRNA_Bc_mockAdj_strictNoTo <- sRNA_Bc_mockAdj_strict %>% anti_join(TomNoBcInBcOnly,by="SEQ")
write_tsv(sRNA_Bc_mockAdj_strictNoTo,"results/Botrytis_matchAligned.mockToFilter_max5ToMock_min5BcOnly.tsv.gz")
