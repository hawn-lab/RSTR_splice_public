library(tidyverse)
library(ggrepel)
library(ggupset)
library(patchwork)
library(cowplot)

#### Volcano plots (A,B) ####
##### Data #####
sleuth_results_tx<-
  read.csv("results/sleuth/SO_adv60_PC_fourCondition_sleuth_table_tx_full_results.csv",
           row.names = "X")%>%
  filter(model_term=="condition" & effect=="RSTR")%>%
  mutate(p_group = case_when(qval<0.001 ~ "FDR < 0.001",
                             qval<0.01 & qval >=0.001 ~ "FDR < 0.01",
                             qval<0.05 & qval >=0.01 ~ "FDR < 0.05", 
                             qval<0.2 & qval>=0.05 ~"FDR < 0.20", 
                             qval<0.3 & qval>=0.2 ~"FDR < 0.30",
                             TRUE ~ "non significant"))%>%
  arrange(-qval)

##### Add gene labels ####
## Find top 25 up/down for labeling in TB condition
tb_top <- sleuth_results_tx %>% 
  filter(reference_level == "LTBI_TB" & b > 0) %>% 
  slice_min(qval, n=25) %>% 
  mutate(label = ext_gene)
tb_bot <- sleuth_results_tx %>% 
  filter(reference_level == "LTBI_TB" & b < 0) %>% 
  slice_min(qval, n=25) %>% 
  mutate(label = ext_gene)
#Also label MEDIA, just do qval < 0.3 b/c there are few enough
media_top <- sleuth_results_tx %>% 
  filter(reference_level == "LTBI_MEDIA" & qval < 0.3) %>% 
  mutate(label = ext_gene)

## Add label info to df
sleuth_results_tx_label<-
  sleuth_results_tx%>%
  full_join(bind_rows(tb_top,tb_bot, media_top)) %>% 
  mutate(label_group = case_when(b < 0 ~ "down",
                                 b > 0 ~ "up")) %>% 
  mutate(facet_label =  case_when(reference_level == "LTBI_MEDIA" ~ "MEDIA \nLTBI <---    ---> RSTR",
                                  reference_level == "LTBI_TB" ~ "Mtb\nLTBI <---    ---> RSTR"))

##### Volcano plot #####
volcano_tx <- sleuth_results_tx_label %>% 
  arrange(-qval) %>% 
  
  ggplot(aes(y=log2(mean_obs), x=b))+
  geom_point(alpha=1, aes(color=p_group))+
  scale_color_manual(values=c("non significant"="grey",
                              "FDR < 0.30" = "black",
                              "FDR < 0.20" = "yellow",
                              "FDR < 0.05" = "darkorange1",
                              "FDR < 0.01" = "red", 
                              "FDR < 0.001"="darkred"))+
  
  # Add repelled labels for down regulated
  geom_text_repel(data=filter(sleuth_results_tx_label,
                              !is.na(label) & label_group == "down"), 
                  aes(label = label), size = 3,
                  nudge_x = -2.2 - filter(sleuth_results_tx_label, 
                                        !is.na(label) & label_group == "down")$b, 
                  segment.size = 0.2, angle = 0, hjust = 0, box.padding = 0.5,
                  segment.color = "grey50", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32) +
  # Add repelled labels for up regulated
  geom_text_repel(data=filter(sleuth_results_tx_label, 
                              !is.na(label) & label_group == "up"), 
                  aes(label = label), size = 3,
                  nudge_x = 2 + filter(sleuth_results_tx_label, 
                                       !is.na(label) & label_group == "up")$b, 
                  segment.size = 0.2, angle = 0, hjust = 0, box.padding = 0.5,
                  segment.color = "grey50", direction = "y",
                  show.legend = FALSE, max.overlaps = Inf, seed=32) +
  geom_vline(xintercept = 0)+
  labs(y="Mean Log2( TPM )", x="Bias Adj Fold Change",
       title = "A)                                                                                                B)")+
  theme_bw()+ 
  guides(color = guide_legend(override.aes = list(alpha=1)))+
  theme(aspect.ratio = 1,
        legend.title = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(vjust = -8)) +
  facet_wrap(~facet_label) +
  lims(x=c(-2.1,2.1))
# volcano_tx

#### Upset plot ####
##### Data #####
DET <- 
  read.csv("results/sleuth/sleuth_model_results/four_condition_contrasts/wald_tests/sleuth_mod_fourCondition_sleuth_table_tx_full_results.csv", 
           row.names = "X") %>% 
  filter(qval < 0.05 & effect=="RSTR")

bulkDEG_interaction<-
  read.csv("results/comp_to_bulkRNASeq/RSTR.interaction_age.sex.batch.model.results.csv.gz")%>%
  filter(variable=="conditionTB:Sample_GroupRSTR"& FDR<0.2)%>%
  select(gene) %>% 
  rename(ext_gene=gene) %>% 
  mutate(reference_level = "lmekin")

##### Signif genes #####
DET.DEG <- DET %>% 
  select(reference_level, b, ext_gene) %>% 
  mutate(b_group = case_when(b < 0 ~ "LTBI",
                             b > 0 ~ "RSTR")) %>% 
  bind_rows(bulkDEG_interaction) %>% #Pretty labels
  mutate(label = case_when(reference_level == "LTBI_MEDIA" & b_group == "LTBI" ~
                             "Sleuth DET genes (media): higher in LTBI",
                           reference_level == "LTBI_MEDIA" & b_group == "RSTR" ~
                             "Sleuth DET genes (media): higher in RSTR",
                           reference_level == "LTBI_TB" & b_group == "LTBI" ~
                             "Sleuth DET genes (TB): higher LTBI expression",
                           reference_level == "LTBI_TB" & b_group == "RSTR" ~
                             "Sleuth DET genes (TB): higher RSTR expression",
                           reference_level == "lmekin" ~
                             "lmekin interaction DEGs"))

gene_df <- DET.DEG %>% 
  #collapse into list column
  group_by(ext_gene) %>% 
  summarise(labels=list(unique(label)), .groups = "drop") 
  
##### plot #####
upset_tx<-
  gene_df %>%
  ggplot(aes(x=labels)) +
  geom_bar() +
  scale_x_upset() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  theme_classic() +
  labs(x="", y="Interection size", title = "C)") + 
  theme(plot.margin = unit(c(1,1,1,175), "pt"))

# Find and save overlaps
gene_wide<-gene_df %>% 
  rownames_to_column() %>% 
  pivot_wider(names_from = labels, values_from = ext_gene) %>% 
  column_to_rownames()

## DET-TB-LTBI only (38)
tb_ltbi <- DET.DEG$ext_gene[!DET.DEG$ext_gene %in% 
                              filter(DET.DEG, label != "Sleuth DET genes (TB): higher LTBI expression")$ext_gene] %>% 
  unique() %>% sort()
length(tb_ltbi) #Check result

## DET-TB-RSTR only (13)
tb_rstr <- DET.DEG$ext_gene[!DET.DEG$ext_gene %in% 
                              filter(DET.DEG, label != "Sleuth DET genes (TB): higher RSTR expression")$ext_gene] %>% 
  unique() %>% sort()
length(tb_rstr) #Check result

## DET-TB-RSTR & DEG lmekin (10)
temp1 <- DET.DEG %>% 
  filter(label == "Sleuth DET genes (TB): higher RSTR expression") %>% 
  pull(ext_gene) %>% unique()
temp2 <- DET.DEG %>% 
  filter(label == "lmekin interaction DEGs") %>% 
  pull(ext_gene) %>% unique()
tb_rstr_kin <- sort(intersect(temp1,temp2))
length(tb_rstr_kin) #Check result

## DET-TB-LTBI & DEG lmekin (7)
temp1 <- DET.DEG %>% 
  filter(label == "Sleuth DET genes (TB): higher LTBI expression") %>% 
  pull(ext_gene) %>% unique()
temp2 <- DET.DEG %>% 
  filter(label == "lmekin interaction DEGs") %>% 
  pull(ext_gene) %>% unique()
tb_ltbi_kin <- sort(intersect(temp1,temp2))
length(tb_ltbi_kin) #Check result

##Format to data frame
list(
  "DET_TB_LTBI_only" = tb_ltbi,
  "DET_TB_RSTR_only" = tb_rstr,
  "DET_TB_LTBI_and_lmekin" = tb_ltbi_kin,
  "DET_TB_RSTR_and_lmekin" = tb_rstr_kin
) %>% 
  plyr::ldply(rbind) %>% 
  column_to_rownames(".id") %>% 
  t() %>% as.data.frame() %>% 
  write_csv("publication/upset_gene_lists.csv", na = "")
## Add plot with total set sizes like UpSetR
# side_plot <- gene_df %>%
#   select(labels) %>%
#   unnest(cols = labels) %>%
#   count(labels) %>%
#   mutate(labels = fct_reorder(as.factor(labels), n)) %>%
#   ggplot(aes(y = n, x = labels)) +
#   geom_col() +
#   coord_flip() +
#   scale_y_reverse() +
#   labs(y = "Set size") +
#   theme_classic() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y = element_blank()) + 
#   theme(plot.margin = unit(c(0,155,-2,0), "pt"))
# 
# upset_tx_all <- cowplot::plot_grid(
#   cowplot::plot_grid(NULL, side_plot, ncol = 1, rel_heights = c(9, 1)),
#   upset_tx, nrow = 1, rel_widths = c(2, 4))
# upset_tx_all

#### Save ####
ggsave("publication/Fig2.DET.volcano.png", volcano_tx,
       width=12, height = 6)
ggsave("publication/Fig2.DET.volcano.tiff", volcano_tx,
       width=12, height = 6)
ggsave("publication/Fig2.DET.upset.png", upset_tx,
       width=5, height = 6)
ggsave("publication/Fig2.DET.upset.tiff", upset_tx,
       width=5, height = 6)
