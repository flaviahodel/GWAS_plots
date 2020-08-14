# Libs
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggthemes)
library(ggrepel)
library(dplyr) # do not load plyr lib afterwards


chicago_plot <- function(df1, df2, threshold) {
  
  message(sprintf("%d SNPs from df1", nrow(df1)))
  
  ### Load dataframe 1 ###
  plot_data_1 <- df1 %>%   
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=as.numeric(max(BP))) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df1, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot))
  
  message(sprintf("%d SNPs from df2", nrow(df2)))
  
  ### Load dataframe 1 ###
  plot_data_2 <- df2 %>%   
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=as.numeric(max(BP))) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df2, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot))
  
  ### -log10 and +log10-transform P-values ###
  plot_data_1$P <- -log10(plot_data_1$P)
  plot_data_2$P <- log10(plot_data_2$P)
  
  ### Genome-wide significant top SNPs per chromosome ###
  tmp_1 <- plot_data_1[plot_data_1$P > -log10(threshold), ]
  lead_snp_1 <- tmp_1 %>% 
    group_by(CHR) %>% 
    slice(which.max(P))
  
  tmp_2 <- plot_data_2[plot_data_2$P < +log10(threshold), ]
  lead_snp_2 <- tmp_2 %>% 
    group_by(CHR) %>% 
    slice(which.min(P))
  
  ### Bind ###
  plot_data <- rbind(plot_data_1, plot_data_2)
  
  ### Generate x-axis ###
  axisdf <- plot_data %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum)) / 2 )
  
  ### Create plot ###
  plot <- ggplot(plot_data, aes(x=BPcum, y=P)) + 
    #specify the y and x values
    geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 1.3) + 
    # create scatterplot colored by chromosome
    scale_color_manual(values = rep(c("#4E79A7", "#A0CBE8"), 22)) + 
    # set a colour pattern 
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center, expand = c(0.03, 0.03)) + 
    # scale the x-axis
    scale_y_continuous(limits = c(-50,50), expand = expansion(mult = 0, add = 0)) + 
    # add x label
    xlab("Chromosome") + 
    # add y label
    ylab(expression(+log[10](italic(P))~'\t'~'\t'~'\t'~-log[10](italic(P)))) + 
    # add annotation value
    geom_label_repel(data = lead_snp_1, 
                     aes(label = SNP), 
                     size = 4,
                     nudge_y = 2,
                     segment.size = 0.4,
                     vjust = -1.5,
                     hjust = 1.5) + 
    # add annotation value
    geom_label_repel(data = lead_snp_2, 
                     aes(label = SNP), 
                     size = 4,
                     nudge_y = 2,
                     segment.size = 0.4,
                     vjust = 1.5,
                     hjust = 1.5) +
    geom_point(data = lead_snp_1, color = "red", size = 2) + # Add highlighted points 
    geom_point(data = lead_snp_2, color = "red", size = 2) + # Add highlighted points 
    geom_hline(yintercept = c(0, -log10(threshold), log10(threshold)), linetype = c("solid", "dashed", "dashed")) + 
    theme_light() +
    theme(legend.position="none",
          plot.margin = margin(6,6,-2,6),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          panel.background  = element_blank())
  
  ### Return the final plot ###
  return(plot)
}
