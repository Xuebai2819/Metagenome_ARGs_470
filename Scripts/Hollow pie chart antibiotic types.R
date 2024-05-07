################### Hollow pie chart ###############
#
library(dplyr)
# Grouped sums using dplyr
for (i in c("chicken", "human", "pig", "soil")) {# i="chicken"
  data <- read.csv(paste0("Antibiotic_class_",i,"_identity04_counts.csv"), header = T)
  result <- data %>%
    group_by(Antibiotic_class) %>%
    summarize(across(everything(), sum, na.rm = TRUE))
  write.csv(result, paste0("Antibiotic_class_",i,"_identity04_counts_unique.csv"),row.names = FALSE)
}
#
for (i in c("chicken", "human", "pig", "soil")) {
  all.data <- read.csv(paste0("Antibiotic_class_",i,"_identity04_counts_unique.csv"), row.names = 1, header = T)
  #
  df <- t(all.data)
  #
  data.prop <- df/rowSums(df)
  #
  relative <- t(data.prop)
  #
  relative_ave <- cbind(relative, rowave=rowMeans(relative))
  #
  relative_ave <- as.data.frame(relative_ave)
  #
  order_relative_ave <- relative_ave[order(relative_ave$rowave,decreasing = T),]
  write.csv(order_relative_ave, paste0("Antibiotic_class_",i,"_identity04_relative_unique.csv"))
}
rm(list=ls())

# Count the number of repeated Antibiotic_class names and merge them with the relative abundance values
#
for (i in c("chicken", "human", "pig", "soil")) {# i="chicken"
  df1 <- read.csv(paste0("Antibiotic_class_",i,"_identity04_counts.csv"), header = T)
  df2 <- read.csv(paste0("Antibiotic_class_",i,"_identity04_relative_unique.csv"), header = T)
  #
  colnames(df2)[1] <- "Antibiotic_class"
  #
  df2_1 <- df2[,c("Antibiotic_class","rowave"), drop = FALSE]
  # Calculate number
  # Sum Antibiotic_class across all samples, i.e. sum the rows
  df1_rowsum <- cbind(df1,rowsum = rowSums(df1[, -1]))
  # 
  df1_rowsum_0 <- df1_rowsum[df1_rowsum$rowsum != 0, ]
  # Use the table function to count the number of repeated elements
  result_df <- data.frame(table(df1_rowsum_0$Antibiotic_class))
  # 
  colnames(result_df) <- c("Antibiotic_class", "Numbers")
  # 
  merged_df <- merge(df2_1, result_df, by = "Antibiotic_class", all = TRUE)
  #
  order_merged_df <- merged_df[order(merged_df$rowave,decreasing = T),]
  write.csv(order_merged_df, paste0("Antibiotic_class_",i,"_identity04_relative_numbers_unique.csv"), row.names = FALSE)
}
rm(list=ls())

# plot
library(ggplot2)
library(patchwork)
library(ggrepel)
library(stringr)
#
pie_top15_list_Abundance <- list()
pie_top15_list_Number <- list()
mycols <- c("#358180", "#3C468E", "#DB381B", "#5B1776", '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A',"#ADB6B699")
for (i in c("chicken", "human", "pig", "soil")) {# i="chicken"，i ="human", i ="pig", i="soil"
  df <- read.csv(paste0("Antibiotic_class_",i,"_identity04_relative_numbers_unique.csv"), row.names = 1, header = T)
  # 
  df_subset <- df[1:15, ]
  # 
  sum_of_remaining_rows <- colSums(df[16:nrow(df), ])
  # 
  sum_of_remaining_rows_df <- t(as.data.frame(sum_of_remaining_rows))
  # 
  rownames(sum_of_remaining_rows_df)[1] <- "Others"
  # 
  merged_df <- rbind(df_subset, sum_of_remaining_rows_df)
  merged_df$Antibiotic_class <- rownames(merged_df)
  # 
  merged_df$rowave <- round(merged_df$rowave, 3)
  # 
  colnames(merged_df) <- c("Abundance", "Number", "Antibiotic_class")
  # 
  merged_df[["Antibiotic_class"]] = factor(merged_df[["Antibiotic_class"]], levels = merged_df[["Antibiotic_class"]], ordered = TRUE)
  #
  merged_df_sum <- merged_df %>%
    mutate(cumsum_Abun = cumsum(Abundance)) %>%
    mutate(cumsum_Num = cumsum(Number))
  #
  merged_df_sum_percent <- merged_df_sum %>%
    mutate(Percentage_num = round(rowSums(select(., starts_with("Number")))/sum(rowSums(select(., starts_with("Number")))), 3) * 100)
  #Calculate the cumulative percentage of number
  merged_df_sumpercent <- merged_df_sum_percent %>%
    mutate(sumpercent_num = cumsum(Percentage_num))
  #
  write.csv(merged_df_sum_percent, paste0("pie_ready_",i,".csv"))
  # Loop to save images into a list
  pie_top15_list_Number[[i]] <- ggplot(data = merged_df_sumpercent, aes(x = 1, y = Percentage_num)) +
    geom_col(color = "white", width = 1, aes(fill = Antibiotic_class)) +
    scale_fill_manual(values = mycols)+
    coord_polar(theta = "y")+
    xlim(-0.5, 1.5)+
    geom_text(data = filter(merged_df_sumpercent, Percentage_num >= 7),
              size = 4,
              color = "white",
              aes(y = 100 - sumpercent_num + 0.5 * Percentage_num,
                  label = str_c(Percentage_num, "%"))) +
    labs(x = "", y = "", title = i) +   
    theme_test() + 
    theme(axis.text = element_blank()) + 
    theme(axis.ticks = element_blank()) +   
    theme(panel.grid=element_blank()) +    
    theme(panel.border=element_blank()) +  
    theme(legend.title = element_blank(), 
          legend.position = "right", 
          axis.title = element_text(size = 12), 
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(5, "pt"),
          plot.title = element_text(size = 12,face = "bold",family = 'serif'))  
  #
  pie_top15_list_Abundance[[i]] <- ggplot(data = merged_df_sum, aes(x = 1, y = Abundance)) +
    geom_col(color = "white", width = 1, aes(fill = Antibiotic_class)) +
    scale_fill_manual(values = mycols)+
    coord_polar(theta = "y")+
    xlim(-0.5, 1.5)+
    geom_text(data =  filter(merged_df_sum, Abundance >= 0.03),
              size = 4,
              color = "white",
              aes(y = 1 - cumsum_Abun + 0.5 * Abundance,
                  label = str_c(Abundance*100, "%"))) +
    labs(x = "", y = "", title = i) +   
    theme_test() + 
    theme(axis.text = element_blank()) +   
    theme(axis.ticks = element_blank()) +   
    theme(panel.grid=element_blank()) +    
    theme(panel.border=element_blank()) +   
    theme(legend.title = element_blank(), 
          legend.position = "right", 
          axis.title = element_text(size = 12), 
          legend.key.height = unit(5, "pt"),
          legend.key.width = unit(5, "pt"),
          plot.title = element_text(size = 12,face = "bold",family = 'serif'))  
}
# CombinePlots(plist, ncol = 4)
p1 = wrap_plots(pie_top15_list_Abundance) + plot_layout(ncol = 2)
p2 = wrap_plots(pie_top15_list_Number) + plot_layout(ncol = 2)
ggsave("pie_Antibiotic_class_Abundance_top15_all.pdf", plot=p1, width = 10, height = 8)
ggsave("pie_Antibiotic_class_Number_top15_all.pdf", plot=p2, width = 10, height = 8)

# #Clear data and canvas
dev.off()
rm(list=ls())

