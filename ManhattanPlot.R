library(tidyverse) # tidyverse packages
library(RColorBrewer) # complement to ggplot
library(ggrepel) # complement to ggplot
library(kableExtra) # table layout

manh_plot <- function(df, threshold) {
  
  ### 1. Compute the cumulative position of SNP ### 
  plot_data <- df %>%   
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=as.numeric(max(BP))) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=as.numeric(BP+tot))
  
  ### 2. Generate x-axis ###
  axisdf <- plot_data %>% 
    group_by(CHR) %>% 
    summarize(center=(max(BPcum) + min(BPcum)) / 2 )
  
  ### 3. Top SNP  ###
  lead_snp_df <- plot_data[order(plot_data$P),]
  lead_snp <- lead_snp_df[1,]
  
  ### 4. Create plot ###
  plot <- ggplot(plot_data, aes(x=BPcum, y=-log10(P))) + 
    #specify the y and x values
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) + 
    # create scatterplot colored by chromosome
    scale_color_manual(values = rep(c("#4E79A7", "#A0CBE8"), 22)) + 
    # set a colour pattern 
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0.03, 0.03)) + 
    # scale the x-axis
    scale_y_continuous(expand = c(0.01, 0.03)) + 
    # remove space between plot area and x axis
    ylim(0,20) +
    theme_light() +
    theme(legend.position="none",
          plot.margin = margin(6,6,-2,6),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line = element_line(color = "black")) +
    xlab("Chromosome") + 
    # add x label
    geom_label_repel(data=lead_snp, # add annotation value
                      aes(label=SNP), 
                     size=4,
                     nudge_y = 2,
                     segment.size = 0.5,
                     vjust = 2,
                     hjust = 1.5) + # add annotation
    geom_point(data= lead_snp, # add annotation value
               color="red", size=2) + # Add highlighted points 
    geom_hline(yintercept = -log10(threshold), linetype="dashed") # threshold line
  
  return(plot) # return the final plot
}


# Theme publication 
theme_publication <- function(base_size = 16, base_family = "") {
  requireNamespace('ggplot2')
  thm = theme_bw(base_size = base_size, base_family = base_family)
  thm + theme(axis.line = element_line(),
              panel.border = element_blank(),
              panel.background  = element_blank(),
              panel.grid.major  = element_blank(),
              panel.grid.minor  = element_blank())}
