######################## predict ####################################
##############Core resistome#############
card <- read.csv(paste0("card_ARO_profile_RPKM_0.4.csv"), row.names = 1, header = T)
#Extract AROID and AROname
AROID_name <- data.frame(ARO_name = card[,1], row.names = row.names(card))
#Remove AROname column
card_all <- card[,-1]
#
t_card_all <- t(card_all)
# 
group_table <- read.csv("metadata.csv",row.names = 1, header = T)
# 
t_card_all = t_card_all[rownames(group_table), ]
# 
card_all_group <- cbind(group_table$species, t_card_all)
# 
colnames(card_all_group)[1] <- "species"
card_all_group <- data.frame(card_all_group)
###########################################
# merge
for (i in c("chicken", "human", "pig", "soil")) { # i = "human"
  species_card_all <- card_all_group[card_all_group$species == i, ]
  # 
  species_card_all <- species_card_all[, -1]
  # Write files corresponding to card_tax of different species
  write.csv(species_card_all, paste0("card_ARO_profile_RPKM_0.4_",i,".csv"))
  card_i <- read.csv(paste0("card_ARO_profile_RPKM_0.4_",i,".csv"),row.names = 1, header = T)
  # Sum the abundances of each sample
  card_abun <- data.frame(rowSums(card_i))
  # 
  colnames(card_abun) <- "Abun"
  ## Count the number of drug-resistant genes
  # ill in the expressed value with 1
  species_card_all[species_card_all > 0] = 1
  species_card_all <- data.frame(species_card_all)
  write.csv(species_card_all, paste0("card_ARO_profile_RPKM_0.4_1_",i,".csv"))
  species_ARGs1 <- read.csv(paste0("card_ARO_profile_RPKM_0.4_1_",i,".csv"),row.names = 1, header = T)
  # Calculate the sum of each row
  table_sample_ARGs_row = data.frame(rowSums(species_ARGs1))
  # 
  colnames(table_sample_ARGs_row) <- "ARO_Freq"
  ARO_level <- cbind(card_abun, table_sample_ARGs_row)
  # Save data table
  write.csv(ARO_level, paste0("ARO_level_",i,".csv"))
  #
  card_i$sampleID <- rownames(card_i)
  ARO_level$sampleID <- rownames(ARO_level)
  # 
  merge.data_ARO <- merge(ARO_level,card_i, by = "sampleID")
  #
  data_ARO <- merge.data_ARO[,4:1372] 
  #
  card_i_t <- data.frame(t(card_i[,-1370]))
  # average
  row_means <- data.frame(rowMeans(card_i_t, na.rm = TRUE))
  #
  colnames(row_means) <- "Ave_Abun"
  #
  card_i_t[card_i_t > 0] = 1
  card_i_t <- data.frame(card_i_t)
  write.csv(card_i_t, paste0("card_ARO_profile_RPKM_0.4_1_t_",i,".csv"))
  species_ARGs1_t <- read.csv(paste0("card_ARO_profile_RPKM_0.4_1_t_",i,".csv"),row.names = 1, header = T)
  # Calculate the sum of each row
  table_sample_ARGs_row_t = data.frame(rowSums(card_i_t))
  # 
  colnames(table_sample_ARGs_row_t) <- "Freq"
  ARO_level_t <- cbind(AROID_name, row_means, table_sample_ARGs_row_t)
  # Remove rows with value 0
  ARO_level_t_filtered <- ARO_level_t[ARO_level_t$Ave_Abun != 0, ]
  # Save data table
  write.csv(ARO_level_t_filtered, paste0("ARO_Freq_abun_",i,".csv"))
  #
  gene_name <- read.csv(paste0("ARO_Freq_abun_",i,".csv"), header = TRUE, check.names = FALSE)
  colnames(gene_name)[1] <- "AROID"
  data_ARO <- data_ARO[,match(gene_name$AROID,names(data_ARO))]
  names(data_ARO) <- gene_name$ARO_name
  
  Richness <- ARO_level$ARO_Freq #total ARGs number of each number
  Abun <- rowSums(data_ARO)  #total abundance of each sample
  #
  library(psych)
  #abundance
  single_Corr_Abun <- corr.test(data_ARO,Abun,use = "complete",method = "spearman")
  single_Abun_r <- single_Corr_Abun$r
  single_Abun_p <- single_Corr_Abun$p
  
  #Richness
  single_Corr_Richness <- corr.test(data_ARO,Richness,use = "complete",method = "spearman")
  single_Richness_r <- single_Corr_Richness$r
  single_Richness_p <- single_Corr_Richness$p
  
  #the gene with the highest correction to the total AMR level was selected as initial gene subset.
  single_mean_r <- (single_Abun_r+single_Richness_r)/2
  abudance_richness_table <- data.frame(cbind(single_Abun_r,single_Abun_p,single_Richness_r,single_Richness_p,single_mean_r))
  # 
  abudance_richness_table <- data.frame(ARO_name = rownames(abudance_richness_table), abudance_richness_table)
  # 
  colnames(abudance_richness_table) <- c("ARO_name", "Abun_r", "Abun_p", "Richness_r", "Richness_p", "merge_r")
  write.csv(abudance_richness_table, paste0("abudance_richness_rp_",i,".csv"),row.names = F)
  single_mean_best_r <- max(single_mean_r)
  #
  single_best_num <- which(single_mean_r== single_mean_best_r)
  single_best_name <- rownames(single_mean_r)[single_best_num]
  #
  single_Abun_best_r <- single_Abun_r[single_best_num]
  single_Abun_best_p <- single_Abun_p[single_best_num]
  single_Rich_best_r <- single_Richness_r[single_best_num]
  single_Rich_best_p <- single_Richness_p[single_best_num]
  #
  #the initial gene subset combined with every other gene and calculated the correction between the sum abundance of each pair gene and total AMR level to find the best combination of pair, and the gene subset was updated
  library(psych)
  select <- single_best_num
  final_result <- c(single_best_name,single_Abun_best_r,single_Abun_best_p,single_Rich_best_r,single_Rich_best_p,single_mean_best_r)
  for (j in 1:9) {
    result <- NULL
    seq <- seq(1,ncol(data_ARO),1)
    for (k in seq[-select]) {
      sum <- rowSums(data_ARO[c(select,k)])
      Abun_corr <- corr.test(sum,Abun,use = "complete",method = "spearman")
      Abun_cor <- Abun_corr$r
      Abun_p <- Abun_corr$p
      Rich_corr <- corr.test(sum,Richness,use = "complete",method = "spearman")
      Rich_cor <- Rich_corr$r
      Rich_p <- Rich_corr$p
      merge_r <- (Abun_cor+ Rich_cor)/2
      a <- c(names(data_ARO)[k],Abun_cor,Abun_p,Rich_cor,Rich_p,merge_r)
      result <- rbind(result,a)
    }
    result <- as.data.frame(result,stringsAsFactors = FALSE)
    names(result) <- c("gene","Abun_r","Abun_p","Richness_r","Richness_p","merge_r")
    
    result[,2:6] <- lapply(result[,2:6], as.numeric)
    r.best <- max(result$merge_r)
    r.name <- result$gene[which(result$merge_r== r.best)] 
    num <- which(names(data_ARO)== r.name)
    best.data <- result[which(result$merge_r== r.best),]
    final_result <- rbind(final_result,best.data)
    select <- c(select,num)
  }
  write.csv(final_result, paste0("predict_Abun_Richness_",i,".csv"))
  }
###################################################################################
# plot
library(readxl)
# 
file_path <- "abundance_richness_best.xlsx"
# 
sheets <- excel_sheets(file_path)
# 
data_list <- list()
# Iterate through all worksheet names and read the data of each worksheet
for (sheet in sheets) { # sheet = "human"
  # 
  final_result <- read_excel(file_path, sheet = sheet)
  final_result$r <-round(final_result$r ,digits = 3)
  final_result$group <- factor(final_result$group, levels = c("Abundance","Richness","Combination"))
  final_result$ID <- as.factor(final_result$ID)
  #
  P <- ggplot(final_result, aes(x=final_result$ID, y=final_result$r, group=group,color=group)) +
    geom_line(linetype="dotted") +
    #scale_color_manual(values = c("#D87070","#3288BD","#469393"))+ 
    geom_point(size=2, shape=20) +
    geom_text(aes(label=gene),nudge_x = 0.15, nudge_y = -0.01,size = 3)+
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.85,0.14),
          legend.box.background = element_rect(color="grey", size=0.5))+
    labs(x="Number of AMR genes",y="Spearman correlation")
  ggsave(paste0("predictor_best_",sheet,".pdf"), plot=P, width =6, height = 5)
}

# Clear data and canvas
dev.off()
rm(list=ls())
