library(tidyverse)
library(BIGpicture)
library(patchwork)

#### Data ####
sleuth_results_tx<-
  read.csv("results/sleuth/sleuth_mod_fourCondition_sleuth_table_tx_full_results.csv",
           row.names = "X")%>%
  filter(model_term=="condition" & effect=="RSTR") 

tx_signif <- sleuth_results_tx %>% 
  filter(qval < 0.05) %>% 
  distinct(ext_gene) %>% 
  pull(ext_gene)

#### Map to string ####
map <- map_string(tx_signif)

#Not mapped
tx_signif[!tx_signif %in% map$map$gene]

#### Plot ####
#Connected genes
## discard parameter depreciated in next update of BIGpicture. Use edge_min and edge_max instead
plot1 <- plot_string(map, node_size = 0.3,
                    discard = "orphan")

#Orphan genes (no connections)
orphans <- tx_signif[!tx_signif %in% plot1$data$symbol]
map2 <- map_string(sort(orphans))
plot2 <- plot_string(map2, node_size = 0.3,
                    layout = "grid")

plot_all <- plot1 / plot2 + plot_layout(heights = c(1.5,1))
# plot_all

#### Save ####
ggsave("publication/FigS2.string.pdf", plot_all, 
       width=6, height=7)
