chicago_plot <- function(df1, df2, threshold) {
  
  requireNamespace('ggplot2')
  message(sprintf("Plotting %d points...", nrow(data)))
  
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
  
  ### Log10
  plot_data_1$P <- -log10(plot_data_1$P)
  plot_data_2$P <- log10(plot_data_2$P)
  
  ### 3. Top SNP  ###
  tmp_1 <- plot_data_1[plot_data_1$P > -log10(threshold), ]
  lead_snp_df_1 <- tmp_1[order(tmp_1$P, decreasing = TRUE), ]
  lead_snp_1 <- lead_snp_df_1[1,]
  
  tmp_2 <- plot_data_2[plot_data_2$P < log10(threshold), ]
  lead_snp_df_2 <- tmp_2[order(tmp_2$P, decreasing = FALSE), ]
  lead_snp_2 <- lead_snp_df_2[1,]
  
  # Bind
  plot_data <- rbind(plot_data_1, plot_data_2)
  
  ### 2. Generate x-axis ###
  axisdf <- plot_data %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum)) / 2 )
  
  ### 4. Create plot ###
  ggplot(plot_data, aes(x=BPcum, y=P)) + 
    #specify the y and x values
    geom_point( aes(color=as.factor(CHR)), alpha = 0.8, size = 1.3) + 
    # create scatterplot colored by chromosome
    scale_color_manual(values = rep(c("#4E79A7", "#A0CBE8"), 22)) + 
    # set a colour pattern 
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center, expand = c(0.03, 0.03)) + 
    # scale the x-axis
    scale_y_continuous(limits = c(-max(abs(plot_data$P)+2), max(abs(plot_data$P)+2)), expand = expansion(mult = 0, add = 0)) + 
    # remove space between plot area and x axis
    theme_publication() + 
    theme_light() +
    theme(legend.position="none",
          axis.text.x = element_text(size = rel(0.75)),
          plot.margin = margin(6,6,-2,6),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line = element_line(color = "black")) +
    # add x label
    xlab("Chromosome") + 
    # add y label
    ylab(expression(+log[10](italic(P))~'\t'~'\t'~'\t'~-log[10](italic(P)))) + 
    # add annotation value
    geom_label_repel(data = lead_snp_1, 
                     aes(label = SNP), 
                     size = 4,
                     nudge_y = 2,
                     segment.size = 0.5,
                     vjust = 2,
                     hjust = 1.5) + 
    # add annotation value
    geom_label_repel(data = lead_snp_2, 
                     aes(label = SNP), 
                     size = 4,
                     nudge_y = 2,
                     segment.size = 0.5,
                     vjust = 2,
                     hjust = 1.5) 
  geom_point(data = lead_snp_1, 
             color = "red", size = 2) + # Add highlighted points 
    geom_point(data = lead_snp_2, # add annotation value
               color = "red", size = 2) + # Add highlighted points 
    geom_hline(yintercept = c(0, -log10(threshold), log10(threshold)), linetype = "dashed") 

  return(plot) # return the final plot
}
