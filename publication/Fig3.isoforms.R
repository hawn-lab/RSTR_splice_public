library(tidyverse)
library(rtracklayer)
library(patchwork)
library(biomaRt)

#### Data ####
#Genes to plot
genes.OI <- c("PDE4A","ZEB2","TIMP1","ACSL4","LST1","GAPDH")
#Significant cutoff of transcripts
sig <- 0.05

#Transcript results
sleuth_results_tx <- read.csv("results/sleuth/SO_adv60_PC_fourCondition_sleuth_table_tx_full_results.csv", 
                              row.names = "X") %>%
  filter(effect == "RSTR")%>%
  mutate(stim = ifelse(grepl("_TB", reference_level), "Mtb",
                       "Media")) %>% 
  filter(ext_gene %in% genes.OI)

#Ref genome
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "uswest")

## Remove decimal from ID
ID.OI <- strsplit(unique(sleuth_results_tx$target_id), split="[.]") %>% 
  sapply(`[`, 1)

tx.key <- getBM(attributes = c("ensembl_transcript_id","external_transcript_name"),
                filters = "ensembl_transcript_id",
                values = ID.OI,
                mart = ensembl)

## Add key names to tx results
sleuth_results_tx_anno <- sleuth_results_tx %>% 
  separate(target_id, into="ensembl_transcript_id", remove = FALSE, extra="drop") %>% 
  full_join(tx.key)

#### plot genes of interest ####
plot.ls <- list()

for(i in genes.OI){
  print(i)
  #### Data ####
  #Subset to gene
  res_tx_sub <- sleuth_results_tx_anno %>% 
    filter(ext_gene == i)
  
  # Assign color code to significant transcripts
  res_tx_sub<-
    res_tx_sub%>%
    group_by(external_transcript_name)%>% # group by Tx
    mutate(tx_color_code = ifelse(min(pval) < sig, external_transcript_name, 
                                  "P > 0.05"))%>%  # code as non-significant if no pval<0.05
    ungroup() %>% # ungroup before ordering color codes
    # Order color coding smartly
    mutate(tx_color_code = factor(tx_color_code,
                                  levels= c("P > 0.05", 
                                            sort(unique(res_tx_sub$external_transcript_name))))) %>% 
    droplevels()
  
  #### Color palette ####
  plot_palette<-c("grey", 
                  '#1f78b4','#33a02c','#ff7f00','#6a3d9a','#e31a1c',
                  '#a6cee3','#b2df8a','#fdbf6f','#cab2d6','#fb9a99') 
  # assign tx names to palette
  all_tx <- levels(res_tx_sub$tx_color_code)
  if("P > 0.05" %in% all_tx){
    names(plot_palette) <- all_tx
  }else{
    names(plot_palette) <- c(NA, all_tx)
  }
  # drop unused colors
  plot_palette<-plot_palette[-which(is.na(names(plot_palette)))] 
  
  #### Plot ####
  tx_plot<-
    res_tx_sub%>%
    arrange(tx_color_code)%>%
    
    ggplot(aes(x=-log10(pval), y=b, group=stim, color=tx_color_code))+
    geom_hline(yintercept = 0, color="grey",
               linetype="dashed", size=0.25)+
    geom_vline(xintercept = -log10(0.05), color="grey",
               linetype="dashed", size=0.25)+
    geom_errorbar(aes(ymin = b-se_b, ymax = b+se_b), width=0, 
                  position=position_dodge(width=0.25))+
    geom_point(aes(size=mean_obs, shape=stim),alpha=0.5,
               position=position_dodge(width=0.25),
               stroke=1)+
    scale_shape_manual(values = c("Media" = 1, "Mtb" = 16))+
    scale_color_manual(values = plot_palette)+
    scale_size(limits = c(0.5,10), breaks = c(1,5,9)) +
    lims(x=c(0,NA)) +
    ylab("Bias Adj Fold Change")+
    xlab("-log10( DET P-value )")+
    labs(shape = " ", 
         color = "Transcript", 
         size = "Log2 Mean\nNormalized\nExpression")+
    ggtitle(i)+
    theme_bw()+
    theme(aspect.ratio=1, 
          panel.grid = element_blank(), 
          strip.background = element_rect(fill="white"),
          plot.title = element_text(hjust=0.5))+
    facet_wrap(~stim)
  
  #Remove some legends
  tx_plot2 <- tx_plot +
    guides(size="none", shape="none")
  plot.ls[[i]] <- tx_plot2
  
  #Save size and shape legends to add back
  if(i == genes.OI[length(genes.OI)]){
    lengend.plot <- tx_plot +
      guides(color = "none",
             shape = guide_legend(override.aes = list(size = 5)),
             size = guide_legend(title.position = "left")) +
      theme(legend.position = "bottom", legend.direction = "vertical")
    
    plot.ls[["space1"]] <- plot_spacer()
    plot.ls[["legend"]] <- cowplot::get_legend(lengend.plot)
    plot.ls[["space2"]] <- plot_spacer()
  }
  
  
}

#### Combine and save ####
plot_all <- wrap_plots(plot.ls) + 
  plot_layout(heights = c(1,1,0.3)) +
  plot_annotation(tag_levels = 
                    list(c("A)","B)","C)","D)","E)","F)"," ")))
# plot_all

ggsave("publication/Fig3.isoform.new.png", plot_all,
       width=15, height=7.5)
ggsave("publication/Fig3.isoform.new.tiff", plot_all,
       width=15, height=7.5)
