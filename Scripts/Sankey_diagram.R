############### Sankey diagram ####################
card_tax_level_freq <- read.csv("card_tax_reads_number.csv", row.names = 1, header = T)
#
card_tax_level_freq_1 <- card_tax_level_freq[,-486]
#
t_card_tax_level_freq_1 <- t(card_tax_level_freq_1)
# 
level_freq_group <- read.csv("metadata.csv",row.names = 1, header = T)
# 
t_card_tax_level_freq_name_1 = t_card_tax_level_freq_1[rownames(level_freq_group), ]
# 
level_freq_df <- cbind(level_freq_group$species, t_card_tax_level_freq_name_1)
# 
colnames(level_freq_df)[1] <- "species"
level_freq_df <- data.frame(level_freq_df)
# 
for (i in c("chicken", "human", "pig", "soil")) {# i="chicken"
  species_level_freq <- data.frame(level_freq_df[level_freq_df$species == i, ])
  # 
  species_level_freq_r <- rbind(level_freq_df[1:16,], species_level_freq)
  # 
  species_level_freq_r <- species_level_freq_r[, -1]
  # 
  species_level_freq_t <- t(species_level_freq_r)
  # 
  write.csv(species_level_freq_t, paste0("card_tax_reads_number_",i,".csv"))
}
species_level_freq_t <- data.frame(species_level_freq_t) 

############# Extract RAGs corresponding to species ##################
card_tax_level_freq_sp <- read.csv("card_tax_reads_number_chicken.csv",row.names = 1, header = T)
#
tax_level_freq <- read.csv("tax.csv", header = T)
# 
index <- tax_level_freq$chicken
#
df <- data.frame()
# 
for (i in index) { #i = "s__Escherichia_coli"
  ff <- card_tax_level_freq_sp[card_tax_level_freq_sp$Species == i,]
  # Remove duplicate rows
  cf <- ff[!duplicated(ff$CARD_ARO_name),]
  df <- rbind(df, cf)
}
write.csv(df, "card_tax_filter_chicken.csv")


card_tax_level_freq_sp <- read.csv("card_tax_reads_number_chicken.csv",row.names = 1, header = T)
tax_level_freq <- read.csv("tax.csv", header = T)
index <- tax_level_freq$chicken
df <- data.frame()
for (i in index) { #i = "s__Escherichia_coli"
  ff <- card_tax_level_freq_sp[card_tax_level_freq_sp$Species == i,]
  #Remove duplicate rows
  cf <- ff[!duplicated(ff$CARD_ARO_name),]
  df <- rbind(df, cf)
}
write.csv(df, "card_tax_filter_chicken.csv")


card_tax_level_freq_sp <- read.csv("card_tax_reads_number_human.csv",row.names = 1, header = T)
tax_level_freq <- read.csv("tax.csv", header = T)
index <- tax_level_freq$human
df <- data.frame()
for (i in index) { #i = "s__Escherichia_coli"
  ff <- card_tax_level_freq_sp[card_tax_level_freq_sp$Species == i,]
  #Remove duplicate rows
  cf <- ff[!duplicated(ff$CARD_ARO_name),]
  df <- rbind(df, cf)
}
write.csv(df, "card_tax_filter_human.csv")


card_tax_level_freq_sp <- read.csv("card_tax_reads_number_pig.csv",row.names = 1, header = T)
tax_level_freq <- read.csv("tax.csv", header = T)
index <- tax_level_freq$pig
df <- data.frame()
for (i in index) { #i = "s__Escherichia_coli"
  ff <- card_tax_level_freq_sp[card_tax_level_freq_sp$Species == i,]
  #Remove duplicate rows
  cf <- ff[!duplicated(ff$CARD_ARO_name),]
  df <- rbind(df, cf)
}
write.csv(df, "card_tax_filter_pig.csv")


card_tax_level_freq_sp <- read.csv("card_tax_reads_number_soil.csv",row.names = 1, header = T)
tax_level_freq <- read.csv("tax.csv", header = T)
index <- tax_level_freq$soil
df <- data.frame()
for (i in index) { #i = "s__Escherichia_coli"
  ff <- card_tax_level_freq_sp[card_tax_level_freq_sp$Species == i,]
  #Remove duplicate rows
  cf <- ff[!duplicated(ff$CARD_ARO_name),]
  df <- rbind(df, cf)
}
write.csv(df, "card_tax_filter_soil.csv")
#
library(data.level_freq)
#install.packages("webshot")
library(webshot)
library(networkD3)
#webshot::install_phantomjs()
# Count the number of ARGs for each species
for (j in c("chicken", "human", "pig", "soil")) { # j = "chicken"
  # 
  filter_level_freq <- read.csv(paste0("card_tax_filter_",j,".csv"),row.names = 1, header = T)
  # Use the table function to calculate the number of repetitions of the species name. The drug-resistant gene name has been de-duplicated, and the number of repetitions is
  result <- data.frame(table(filter_level_freq$Species))
  colnames(result)[1] <- "Species"
  #
  level_name <- filter_level_freq[!duplicated(filter_level_freq$Species),c(2:8)]
  # Merge level and freq tables
  level_freq <- merge(level_name, result, by = "Species")
  level_freq <- level_freq[,c(2:7,1,8)]
  # 
  level_freq_bac <- level_freq[level_freq$Kingdom == "k__unclassified_d__Bacteria", ]
  # Kingdom
  level_freq_bac$Kingdom[level_freq_bac$Kingdom == "k__unclassified_d__Bacteria"] <- "Bacteria"
  #
  level <- colnames(level_freq_bac)[2:7]
  for (l in level) {
    level_freq_bac[, l] <- ifelse(grepl("unclassified", level_freq_bac[, l]), 
                              paste0(l,"_unknown"),
                                gsub(".*__", "", level_freq_bac[, l]))
  }
  # Write a Sankey diagram preparation document
  write.csv(level_freq_bac, paste0("bacteria_level_for_Sankey_",j,".csv"),row.names = F)
}
# plot
#Integrate the nested relation of classification and ARGs number to construct the link list of Sankey graph
# s = "chicken", s = "pig", s = "human", s = "soil"
table_level_freq_bac <- read.csv(paste0("bacteria_level_for_Sankey_",s,".csv"),header = T)
genus_species <- table_level_freq_bac[c('Genus','Species','Freq')]
names(genus_species) <- c('source', 'target', 'Freq')
family_genus <- aggregate(table_level_freq_bac$Freq, by = list(table_level_freq_bac$Family, table_level_freq_bac$Genus), FUN = sum)
names(family_genus) <- c('source', 'target', 'Freq')
order_family <- aggregate(table_level_freq_bac$Freq, by = list(table_level_freq_bac$Order, table_level_freq_bac$Family), FUN = sum)
names(order_family) <- c('source', 'target', 'Freq')
class_order <- aggregate(table_level_freq_bac$Freq, by = list(table_level_freq_bac$Class, table_level_freq_bac$Order), FUN = sum)
names(class_order) <- c('source', 'target', 'Freq')
phylum_class <- aggregate(table_level_freq_bac$Freq, by = list(table_level_freq_bac$Phylum, table_level_freq_bac$Class), FUN = sum)
names(phylum_class) <- c('source', 'target', 'Freq')
kingdom_phylum <- aggregate(table_level_freq_bac$Freq, by = list(table_level_freq_bac$Kingdom, table_level_freq_bac$Phylum), FUN = sum)
names(kingdom_phylum) <- c('source', 'target', 'Freq')

link_list <- rbind(kingdom_phylum,phylum_class, class_order, order_family, family_genus,genus_species)

#Build the node list and assign ID  to the category names in the Link list
node_list <- reshape2::melt(table_level_freq_bac, id = 'Freq')
node_list <- node_list[!duplicated(node_list$value), ]
head(node_list)

link_list$IDsource <- match(link_list$source, node_list$value) - 1 
link_list$IDtarget <- match(link_list$target, node_list$value) - 1
head(link_list)

#plot sankey graph by networkD3 package
#library(networkD3)
color <- 'd3.scaleOrdinal() .domain(["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]) 
.range(["#377EB8", "#8ECFC9","#008080", "#f6a34a", "#ffc482", "#008000","#9D1B45"])'
p <- sankeyNetwork(Links = link_list, Nodes = node_list,
                   Source = 'IDsource', Target = 'IDtarget', Value = 'Freq', 
                   NodeID = 'value', NodeGroup = 'variable',
                   nodeWidth = 30, colourScale = color,
                   fontSize = 10, sinksRight = FALSE,width = 1000,height = 1100)

p
#save as .html
saveNetwork(p,paste0("Bacteria_level_Sankey_",s,".html"))
#library(webshot)
webshot(paste0("Bacteria_level_Sankey_",s,".html"),paste0("Bacteria_level_Sankey_",s,".pdf"),vwidth = 1000,vheight=1100)
# Clear data and close plotting devices
rm(list = ls())
dev.off()
