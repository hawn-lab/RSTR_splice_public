library(tidyverse)

#### Table 1: DET in Media or TB ####
read_csv("results/sleuth/sleuth_mod_fourCondition_sleuth_table_tx_full_results.csv") %>% 
  filter(effect=="RSTR")%>%
  unite("contrast", c(group_2, reference_level), sep="_vs_")%>%
  group_by(contrast)%>%
  summarise(FDR.05 = length(which(qval<0.05)), 
            FDR.1 = length(which(qval<0.1)),
            FDR.2 = length(which(qval<0.2)), 
            FDR.3 = length(which(qval<0.3))) %>% 
  write_csv("publication/Table1.csv")


#### Table 2: Overlap genes with DET and lmekin DEG ####
t2 <- read_csv("results/comp_to_bulkRNASeq/DEG_contrast_compare.csv") %>% 
  group_by(contrast)%>%
  summarise(FDR.05 = length(which(FDR<0.05)), 
            FDR.1 = length(which(FDR<0.1)),
            FDR.2 = length(which(FDR<0.2)), 
            FDR.3 = length(which(FDR<0.3)))

#Nonredundant totals
t2_total <- read_csv("results/comp_to_bulkRNASeq/DEG_contrast_compare.csv") %>% 
  group_by(gene) %>% 
  summarise(min_FDR = min(FDR)) %>% 
  ungroup() %>% 
  summarise(FDR.05 = length(which(min_FDR<0.05)), 
            FDR.1 = length(which(min_FDR<0.1)),
            FDR.2 = length(which(min_FDR<0.2)), 
            FDR.3 = length(which(min_FDR<0.3)))

bind_rows(t2, t2_total) %>% 
  write_csv("publication/Table2.csv")
