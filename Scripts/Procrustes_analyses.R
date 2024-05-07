##########Procrustes analyses###################
#ARGs profile
AMR <- data.frame(t(read.csv("species_filter_ARO_Procrustes_analyses.csv",header=T,row.names = 1)))
AMR$SampleID <- rownames(AMR)
metadata_group <- read.csv("metadata.csv", header = T)
AMR.gene <- merge(metadata_group, AMR, by = "SampleID")
rownames(AMR.gene) <- AMR.gene$SampleID
AMR.gene <- AMR.gene[, -1]  
# Use sepcies as factors
AMR.gene$species <- factor(AMR.gene$species,levels = c("chicken","human","pig","soil"))
# Sort by species
AMR.gene <- AMR.gene[order(AMR.gene$species),]
AMR.gene.Abun <- AMR.gene[,-c(1:2)]

#species data
species.data <- data.frame(t(read.csv("tax_s_all_relative_abundance_0.97_filter_species.csv",header = T,row.names=1))) 

library(vegan)
#Hellinger transform
AMR.gene.Abun.hel <- decostand(AMR.gene.Abun, "hellinger")
species.data.hel <- decostand(species.data, "hellinger")
#
gene.distance <- vegdist(AMR.gene.Abun.hel,method = "bray") 
species.distance <- vegdist(species.data.hel,method = "bray")

# make pcoas 
library(ape)
pcoa_gene <- as.data.frame(pcoa(gene.distance)$vectors)
pcoa_species <- as.data.frame(pcoa(species.distance)$vectors)

# procrustes
# Compares the similarity of two configurations (ordination) and computes a Procrustes match between them
pro <- procrustes(pcoa_gene, pcoa_species)
pro
# Perform a hypothesis test of Procrustes analysis to check whether the similarity between two configurations (ordination) is significant
pro_test <- protest(pcoa_gene, pcoa_species, perm = 999)  
M2 <- round(pro_test$ss,3)  # Deviation sum of squares statistic M2
pro_test$
#########
#
eigen <- sqrt(pro$svd$d)
#
percent_var <- signif(eigen/sum(eigen), 4)*100
# Configuration matrix in Procrustes analysis results object pro
beta_pro <- data.frame(pro$X)
# Information about the rotated second configuration in Procrustes analysis results
trans_pro <- data.frame(pro$Yrot)
#
beta_pro$UserName <- rownames(beta_pro)
#
beta_pro$type <- "AMR gene"
#
beta_pro$species <- AMR.gene$species
#
trans_pro$UserName <- rownames(trans_pro)
trans_pro$type <- "Species level"
trans_pro$species <- AMR.gene$species

colnames(trans_pro) <- colnames(beta_pro)
#
r <- signif(pro_test$t0, 3)
#
pval <- signif(pro_test$signif, 1)
#
plot <- rbind(beta_pro, trans_pro)
# Drawing
library(ggplot2)
library(ggsci)
plot$species <- factor(plot$species,levels = c("chicken","human","pig","soil"))
gene_species <- ggplot(plot) +
  geom_point(size = 1, alpha=0.75, aes(x = Axis.1, y = Axis.2, color = species,shape = type)) +
  scale_color_manual(values = c("chicken" ="#ed6e69" , "human" = "#34ABE0", "pig" = "#269D77", "soil" = "#AA5A00")) +   #point:scale_color_npg()��bar:scale_fill_npg()
  theme_bw()  +
  geom_line(aes(x= Axis.1, y=Axis.2, group=UserName,color=species), alpha = 0.6, linewidth = 0.8) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.position = 'right',
        axis.text = element_text(size=8),
        axis.title = element_text(size=9),
        aspect.ratio = 1) +
  guides(color = guide_legend(ncol = 1)) +
  xlab(paste0("PCoA 1 (",percent_var[1],"%)")) +
  ylab(paste0("PCoA 2 (",percent_var[2],"%)")) +
  labs(title=paste0("r = ",r,", p = ",pval,", M2 = ",M2),size=6)
gene_species
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 10, spaceLegend = 1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize,face = "bold",vjust = 0), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

pdf("AMR.species.protest.bray.PCoA.pdf", width = 11, height = 8.5)
print(addSmallLegend(gene_species))
dev.off()

#The default color of ggplot2
library(scales)
#show_col(hue_pal()(9)) 
#cols <- c("#F8766D","#B79F00","#00BA38","#00BFC4","#619CFF","#F564E3")
#cols <- hue_pal()(9)
cols <- c("chicken" ="#ed6e69" , "human" = "#34ABE0", "pig" = "#269D77", "soil" = "#AA5A00")
show_col(cols)
table(beta_pro$species)
#
cols.data <- rep(cols,c(120,119,120,111))

##residuals plot
#
residuals <- residuals(pro_test)
#
pdf("gene.species.protest.bray.residuals.pdf", width = 8, height = 6)
# 
plot(pro_test, kind = 2, col = cols.data)

# 添加图例
legend("topright", 
       inset = 0.02, 
       box.col = "grey", 。
       legend = c("chicken", "human", "pig", "soil"),
       col = cols, lty = 1, cex = 0.6)

dev.off()
