library(tidyverse)
library(openxlsx)

#### Data ####
#Four condition model
tx_results <- read_csv("results/sleuth/sleuth_mod_fourCondition_sleuth_table_tx_full_results.csv") %>%
  #add rank 
  group_by(effect) %>%
  mutate(Rank = min_rank(qval)) %>% 
  ungroup()

#Interaction model
interact_results <- read_csv("results/sleuth/sleuth_mod_interaction_sleuth_table_tx_full_results.csv")

#Bulk analysis
orig_degs <- read_csv("data_raw/bulk_analysis/RSTR.interaction_age.sex.batch.model.results.csv.gz") %>%
  filter(variable == "conditionTB:Sample_GroupRSTR") %>% 
  #Add DEG designation
  mutate(DEG = ifelse(FDR<0.2, "Y", "N")) %>% 
  select(gene, DEG) %>% 
  rename(ext_gene=gene)
 
#### Table S1A: Media DET ####
A <- tx_results %>% 
  filter(effect == "RSTR" & reference_level == "LTBI_MEDIA" & qval < 0.2) %>% 
  mutate(condition = "MEDIA") %>% 
  mutate(direction = ifelse(b<0, "down","up")) %>% 
  left_join(orig_degs) %>% 
  select(Rank, ext_gene, target_id, condition, direction, qval, DEG) %>% 
  rename(`Gene Name`=ext_gene, `Target ID`=target_id, `Media/Mtb`=condition,
         FDR=qval, `Mtb*RSTR DEG (Y/N)`=DEG) %>% 
  mutate(FDR = signif(FDR, digits = 3))

#### Table S1B: TB DET ####
B <- tx_results %>% 
  filter(effect == "RSTR" & reference_level == "LTBI_TB" & qval < 0.2) %>% 
  mutate(condition = "TB") %>% 
  mutate(direction = ifelse(b<0, "down","up")) %>% 
  left_join(orig_degs) %>% 
  select(Rank, ext_gene, target_id, condition, direction, qval, DEG) %>% 
  rename(`Gene Name`=ext_gene, `Target ID`=target_id, `Media/Mtb`=condition,
         FDR=qval, `Mtb*RSTR DEG (Y/N)`=DEG) %>% 
  mutate(FDR = signif(FDR, digits = 3))


#### Table S1C: Interaction DET ####
C <- interact_results %>% 
  filter(effect == "interaction" & qval < 0.2) %>% 
  mutate(direction = ifelse(b<0, "down","up")) %>% 
  left_join(orig_degs) %>% 
  select(ext_gene, target_id, direction, qval, DEG) %>% 
  rename(`Gene Name`=ext_gene, `Target ID`=target_id, 
         FDR=qval, `Mtb*RSTR DEG (Y/N)`=DEG)

#### Save ####
df_list <- list("S1A"=A,
                "S1B"=B,
                "S1C"=C)

write.xlsx(df_list, 'publication/TableS1.xlsx')
