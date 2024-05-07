library(reshape2)
library(ggplot2)
library(patchwork)
library(tidyverse)
# Read group information
group_table <- read.csv("metadata.csv",row.names = 1, header = T)
for (i in c("chicken","human", "pig", "soil")) { # i = "chicken"
  card_counts <- read.csv(paste0(i,"_tax_relative_p_bac.csv"), row.names = 1, header = T)
  #
  t_card_counts <- t(card_counts)
  # Extract the grouping information of the species corresponding to i
  group_table_s <- group_table[group_table$species == i, ]
  # 
  t_card_counts = t_card_counts[rownames(group_table_s), ]
  # 
  card_counts_group <- cbind(group_table_s$province, t_card_counts)
  # 
  colnames(card_counts_group)[1] <- "province"
  card_counts_group <- data.frame(card_counts_group)
  stacked_area <- list()
  for (j in c("SC", "YN", "GZ", "JS")) { # j = "SC"
    pr_card_counts <- card_counts_group[card_counts_group$province == j, ]
    # 
    pr_card_counts <- pr_card_counts[, -1]
    # 
    pr_card_counts_t <- data.frame(t(pr_card_counts))
    # 
    write.csv(pr_card_counts_t, paste0(i,"_phylum_relative_abundance",j,"_.csv"))
    # 
    pr <- read.csv(paste0(i,"_phylum_relative_abundance",j,"_.csv"), row.names = 1, header = T)
    # average rows
    pr_row_ave <- cbind(pr,ave=rowMeans(pr))
    # Sort by average in descending order
    pr_ave_order <- pr_row_ave[order(pr_row_ave$ave, decreasing = T), ]
    # Extract top5
    pr_top5 <- pr_ave_order[1:5,]
    # 
    pr_6sum <- rbind(pr_top5, Others=rowSums(pr_ave_order[6:nrow(pr_ave_order), ]))
    # 
    pr_6sum_df <- select(pr_6sum, -ave)
    # 
    pr_6sum_df_t <- t(pr_6sum_df)
    #
    pr_6sum_df_1 <- t(pr_6sum_df_t/rowSums(pr_6sum_df_t))
    # 
    write.csv(pr_6sum_df_1, paste0(i,"_phylum_relative_abundance_",j,"_ready.csv"))
    # 
    dat <- read.csv(paste0(i,"_phylum_relative_abundance_",j,"_ready.csv"),header = F, row.names = 1)
    rownames(dat)[1] <- "sampleID"
    new_rownames <- sub(".*_(.*)", "\\1", rownames(dat))
    rownames(dat) <- new_rownames
    dat <- data.frame(t(dat))
    # Add continuous values
    dat$NO <- 1:nrow(dat)
    dat <- dat[ , c(1,7,2:6,8)]
    #按
    dat <- dat[order(dat$Others), ]
    #
    pr_dat <- melt(dat, id = c("NO", "sampleID"))
    pr_dat$value <- as.numeric(as.character(pr_dat$value))
    stacked_area[[j]] <- ggplot(pr_dat, aes(x = NO , y = value, fill = variable)) +
      geom_area() +
      scale_fill_manual(values = c('#B5B5B7', "#EBB86F", '#719160', '#A69184',  '#775F4D', '#82BFB3')) +
      theme(panel.grid = element_blank(), 
            panel.background = element_blank(), 
            axis.ticks = element_line(color = 'black', linewidth = 0.5), 
            axis.text.x = element_blank(), 
            axis.text.y = element_text(color = 'black', size = 9), 
            axis.title.y = element_text(color = 'black', size = 10),
            legend.text = element_text(color = 'black', size = 9), 
            legend.title = element_text(color = 'black', size = 10)) +
      scale_x_continuous( expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = '', y = 'Relative abundance', title = j, fill = 'Phylum')
  }
  # CombinePlots(plist, ncol = 4)
  p1 = wrap_plots(stacked_area) + plot_layout(ncol = 4)
  ggsave(paste0(i,"_phylum_relative_abundance_province.pdf"), plot=p1, width = 15, height = 3)
}

# Clear data and canvas
dev.off()
rm(list=ls())