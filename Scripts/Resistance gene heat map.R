################## heatmap ###################
library(ComplexHeatmap)
library(dplyr)
library(pals)
library(pheatmap) 
library(RColorBrewer) 
library(circlize) 
# 
col.df <- read.csv("metadata.csv",row.names = 1)
#
df <- read.csv("core_resistance_genes_relative.csv",row.names = 1, header = T) %>%
  t() %>% 
  .[rownames(col.df), ] %>%
  scale() %>%
  t()
# Line comment
row.df <- read.csv("Line_comment.csv", row.names = 1)

# plot
col_table = list(
  Resistance_Mechanism = c("antibiotic target alteration/replacement" = "#556B2F", "antibiotic target alteration" = "#B9B9DD", "antibiotic efflux" = "#3288BD", "antibiotic target protection" = "#E4E569", "antibiotic target replacement" = "#827F88", "antibiotic inactivation" = "#B57C82"),
  Drug_Class = c("Multidrug" = "#1F77B4", "Aminocoumarin" = "#FF7F0E", "Tetracycline" = "#1D933A",  "Fosfomycin" = "#87CEEB", "MLS" = "#9467BD", "Fluoroquinolone" = "#8C564B", "Glycopeptide" = "#5b8e7d", "Peptide" = "#7F7F7F", "Phenicol" = "#E4E569", "Pleuromutilin" = "#3288BD", "Mupirocin" = "#FFD700", "Diaminopyrimidine" = "#FFA07A", "Aminoglycoside" = "#98FB98", "Rifamycin" = "#B57C82", "Beta-lactam" = "#134B5F"),
  species = c("chicken" = "#ed6e69" , "human" = "#34ABE0", "pig" = "#269D77", "soil" = "#AA5A00"),
  province = c("SC" = "#DFFBFF", "GZ" = "#F0C9FF", "YN" = "#FFF9DF", "JS" = "#CFFFC9")
)
head(col_table)
# Check the maximum and minimum values
min(df)
max(df)
# Convert to mean
df[df > 4] = 4
pdf("ARGs_top50_heatmap.pdf", width = 10, height = 8)
heatmap(
  df,
  name = "Zscore",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = col.df,
  annotation_row = row.df,
  annotation_colors = col_table,
  col = colorRampPalette(c("#2166AC","#67A9CF","white","#EF8A62","#B2182B"))(100) 
)
# Clear data and canvas
dev.off()
rm(list=ls())
