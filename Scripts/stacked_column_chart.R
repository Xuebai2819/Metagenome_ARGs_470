# Install R package based on CRAN, puzzle, data cleaning, drawing, data filtering
p_list = c("cowplot", "tidyfst", "ggplot2", "dplyr","patchwork")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

# Set color palette
cb_palette <- c("#556B2F", "#B9B9DD", "#3288BD", "#E4E569", "#827F88", "#B57C82", "#ADD8A4", "#5F50A1", "#7ED0DE", "#1D933A", "#CC141D", 
                "#ADC0E3", "#EB9486", "#9BA791", "#134B5F", "#BB7B53", "#EE8354", "#9D1B45", "#C9B014", "#505D75", "#59A691",   '#2F4F4F',
                "#556B2F", "#B9B9DD", "#3288BD", "#E4E569", "#827F88", "#B57C82", "#ADD8A4", "#5F50A1", "#7ED0DE", "#1D933A", "#CC141D",
                "#ADC0E3", "#EB9486", "#9BA791", "#134B5F", "#BB7B53", "#EE8354", "#9D1B45", "#C9B014", "#505D75", "#59A691",   '#2F4F4F')

# Draw the main image
top4_list <- list()
for (j in c("chicken", "pig", "human", "soil")) {
  all <- read.csv(paste0(j,"_top5_pfg_stacked.csv"), header = T)
  # Transform factor implementation sorting
  all$Taxonomy = factor(all$Taxonomy,levels = unique(all$Taxonomy))
  cut_colors <-setNames(cb_palette, levels(all$Taxonomy))
  cut_colors
  # PhylumFirmicutes     PhylumProteobacteria     PhylumActinobacteria       PhylumBacteroidota 
  # "#556B2F"                "#B9B9DD"                "#3288BD"                "#E4E569" 
  # PhylumPhylumOthers   FamilyLactobacillaceae    FamilyEnterococcaceae  FamilyStaphylococcaceae 
  # "#827F88"                "#B57C82"                "#ADD8A4"                "#5F50A1" 
  # FamilyEnterobacteriaceae       FamilyFamilyOthers       GenusLactobacillus   GenusLigilactobacillus 
  # "#7ED0DE"                "#1D933A"                "#CC141D"                "#ADC0E3" 
  # GenusLimosilactobacillus        GenusEnterococcus         GenusGenusOthers                     <NA> 
  #   "#EB9486"                "#9BA791"                "#134B5F"                "#BB7B53" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#EE8354"                "#9D1B45"                "#C9B014"                "#505D75" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#59A691"                "#2F4F4F"                "#556B2F"                "#B9B9DD" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#3288BD"                "#E4E569"                "#827F88"                "#B57C82" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#ADD8A4"                "#5F50A1"                "#7ED0DE"                "#1D933A" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#CC141D"                "#ADC0E3"                "#EB9486"                "#9BA791" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#134B5F"                "#BB7B53"                "#EE8354"                "#9D1B45" 
  # <NA>                     <NA>                     <NA>                     <NA> 
  #   "#C9B014"                "#505D75"                "#59A691"                "#2F4F4F" 
  
  g <- ggplot(data = all, aes(x = factor(level, levels = c("Phylum", "Family", "Genus")), 
                              y = SumAbun, fill = Taxonomy))+
    geom_bar(stat = "identity")+
    xlab("") +ylab("Proportion")+
    labs(title=j)+
    theme(plot.title = element_text(size = 12,face = "bold",family = 'serif',hjust = 0.5))+
    scale_fill_manual(values = cut_colors)+
    theme_bw()+
    # theme(panel.grid = element_blank(),
    #       panel.background = element_rect(colour = "black",
    #                                       linewidth = 1),
    #       axis.ticks = element_line(colour = "black",
    #                                 linewidth = .5,
    #                                 linetype = 1,
    #                                 lineend = 1),
    #axis.text = element_text(colour = "black"))+
    theme(legend.position  =  "none")
  g
  
  #Legneds
  #Phylum
  Phylum <- all[all$level=="Phylum",]
  label =  gsub("Phylum" , "" , Phylum$Taxonomy)
  for (i in 1:length(label)){label[i] = gsub("_[0-9]*", "", label[i])}
  label
  
  
  P_legend <- all %>% filter(level %in% Phylum$level) %>%
    ggplot(aes(level, fill = Taxonomy))+ 
    geom_bar()+
    scale_fill_manual(values  =  cut_colors[names(cut_colors) %in% Phylum$Taxonomy] ,
                      name  =  "Phylum",
                      labels  =  label)+
    
    theme(legend.key.size  =  unit(0.3, 'cm'),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
  P_legend 
  
  #Family
  Family <- all[all$level=="Family",]
  label =  gsub("Family" , "" , Family$Taxonomy)
  for (i in 1:length(label)){label[i] = gsub("_[0-9]*", "", label[i])}
  label
  
  
  F_legend <- all %>% filter(level %in% Family$level) %>%
    ggplot(aes(level, fill = Taxonomy))+ 
    geom_bar()+
    scale_fill_manual(values  =  cut_colors[names(cut_colors) %in% Family$Taxonomy] ,
                      name  =  "Family",
                      labels  =  label)+
    
    theme(legend.key.size  =  unit(0.3, 'cm'),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
  F_legend 
  
  #Genus
  Genus <- all[all$level=="Genus",]
  label =  gsub("Genus" , "" , Genus$Taxonomy)
  for (i in 1:length(label)){label[i] = gsub("_[0-9]*", "", label[i])}
  label
  
  
  G_legend <- all %>% filter(level %in% Genus$level) %>%
    ggplot(aes(level, fill = Taxonomy))+ 
    geom_bar()+
    scale_fill_manual(values  =  cut_colors[names(cut_colors) %in% Genus$Taxonomy] ,
                      name  =  "Genus",
                      labels  =  label)+
    
    theme(legend.key.size  =  unit(0.3, 'cm'),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
  G_legend 
  
  # Combo figures

  p = plot_grid(g,
                plot_grid(get_legend(P_legend),
                          get_legend(F_legend),
                          get_legend(G_legend),
                          ncol  =  1,
                          align  =  "hv"),
                ncol  =  2,align  =  "hv",rel_widths = c(5,3))
  top4_list[[j]] = p

}
  # 
  width = 119
  height = 89
  zoom = 2
  p1 = wrap_plots(top4_list) + plot_layout(ncol = 2)
  ggsave("all_top5_pfg_stacked.pdf", p1, width  =  width * zoom, height  =  height * zoom, units  =  "mm")
  p1
  


# Clear data and canvas
dev.off()
rm(list=ls())
