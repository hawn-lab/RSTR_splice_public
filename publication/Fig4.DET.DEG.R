library(tidyverse)
library(ggrepel)
library(patchwork)

#### Data ####
#lmekin DEG
DEG.contrast <-  read_csv("results/DEG_contrast_compare.csv") %>% 
  mutate(label = recode(label, 
                        "RSTR - LTBI in Mtb"="Mtb\nLTBI <---  ---> RSTR",
                        "RSTR - LTBI in media"="Media\nLTBI <---  ---> RSTR"))
#DET
DET <- 
  read.csv("results/sleuth/sleuth_model_results/four_condition_contrasts/wald_tests/sleuth_mod_fourCondition_sleuth_table_tx_full_results.csv", 
           row.names = "X")

#Format DET data
DET.format <- DET %>% 
  filter(effect == "RSTR" & qval < 0.05) %>% 
  mutate(contrast = case_when(reference_level == "LTBI_TB" ~ "TB_RSTR - TB_LTBI",
                              reference_level == "LTBI_MEDIA" ~ "MEDIA_RSTR - MEDIA_LTBI")) %>% 
  dplyr::select(ext_gene, contrast, qval) %>% 
  rename(FDR=qval, gene=ext_gene) %>% 
  group_by(gene) %>% 
  mutate(DET_no=row_number()) %>%
  group_by(gene, contrast)%>%
  mutate(min_DET_fdr = min(FDR, na.rm = T),
         n_DET = max(DET_no))%>%
  dplyr::select(-c(DET_no, FDR))%>%
  distinct()

# Combine DEG and DET data
DEG.DET <- DEG.contrast%>%
  select(gene, contrast, FDR, estimate) %>% 
  full_join(DET.format, by=c("gene", "contrast")) %>% 
  mutate(label = recode(contrast, 
                        "TB_RSTR - TB_LTBI"="Mtb\nRSTR vs LTBI",
                        "MEDIA_RSTR - MEDIA_LTBI"="Media\nRSTR vs LTBI")) %>% 
  
  filter(!is.na(min_DET_fdr) & !is.na(FDR))%>% 
  #Set gene labels
  mutate(gene_label =ifelse(FDR<0.05 & min_DET_fdr<0.05, gene, NA)) %>% 
  #Create estimate groups
  mutate(est_group = case_when(estimate < -2 ~ "[-3, -2)",
                               between(estimate, -2, -1) ~ "[-2, -1)",
                               between(estimate, -1, -0.5) ~ "[-1, -0.5)",
                               between(estimate, -0.5, 0) ~ "[-0.5, 0)",
                               between(estimate, 0, 0.5) ~ "[0, 0.5)",
                               between(estimate, 0.5, 1) ~ "[0.5, 1)"
  ),
  est_group = factor(est_group, levels=c("[-3, -2)", "[-2, -1)",
                                         "[-1, -0.5)", "[-0.5, 0)",
                                         "[0, 0.5)", "[0.5, 1)"))) %>% 
  mutate(n_DET_fct = factor(n_DET, levels = c(1,2,4,6))) %>% 
  mutate(label_group = case_when(-log10(FDR) > 4 ~ "highY",
                                 -log10(min_DET_fdr) > 1.5 ~ "highX",
                                 TRUE ~ "low")) %>% 
  arrange(min_DET_fdr)

##### Plot #####
## Define colors
col.vec <- c("#08306b","#2171b5","#6baed6","#c6dbef",
             "#fc9272","#cb181d")
#
fdr.plot <- DEG.DET %>% 
  ggplot(aes(y=-log10(FDR), x=-log10(min_DET_fdr), color=est_group))+
  geom_hline(yintercept = -log10(0.05), color="grey", lty="dashed")+
  geom_vline(xintercept = -log10(0.05), color="grey",lty="dashed")+
  annotate("text", x=7, y=-log10(0.05), label = "FDR = 0.05", 
           vjust = -1, size=3) +
  geom_abline(slope = 1, intercept = 0, color="grey")+
  geom_point(aes(size=n_DET_fct))+
  facet_wrap(~label, ncol=1)+
  #Label non-cluttered genes
  ggrepel::geom_text_repel(data=filter(DEG.DET, label_group == "highY"),
                           aes(label=gene_label), size=3,
                           segment.size = 0.2, box.padding = 0.5,
                           show.legend = FALSE, 
                           seed=32)+
  #Label cluttered genes
  ##Left
  geom_text_repel(data=filter(DEG.DET, label_group == "low"),
                  aes(label = gene_label), size = 3,
                  nudge_x = log10(filter(DEG.DET, label_group == "low")$min_DET_fdr),
                  segment.size = 0.2, angle = 0, hjust = 0, 
                  direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32,
                  min.segment.length = 0.1) +
  ##Right
  geom_text_repel(data=filter(DEG.DET, label_group == "highX"),
                  aes(label = gene_label), size = 3,
                  nudge_x = 2.7+log10(filter(DEG.DET, label_group == "highX")$min_DET_fdr),
                  segment.size = 0.2, angle = 0, hjust = 0, 
                  direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32,
                  min.segment.length = 0.1) +
  scale_y_continuous(limits = c(0, NA))+
  scale_x_continuous(limits = c(0, NA))+
  scale_color_manual(values=col.vec) +
  theme_classic() +
  labs(x="-log10( min DET FDR )", y="-log10( lmekin FDR )",
       size="Number of DETs", color="lmekin log2\nfold change") +
  coord_fixed() +
  guides(color = guide_legend(override.aes=list(shape = 15, size=6)))

fdr.plot

#### Save ####
ggsave("publication/Fig4.DET_DEG_alt.png", fdr.plot,
       width=8, height=10)
ggsave("publication/Fig4.DET_DEG_alt.tiff", fdr.plot,
       width=8, height=10)
