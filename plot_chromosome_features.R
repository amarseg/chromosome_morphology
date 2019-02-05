library(tidyverse)
library(ggpubr)
library(gridExtra)

chromosome_data <- read_csv('all_filaments.csv') %>%
  select(-'ID')

variable_name <- colnames(chromosome_data[1:22])

avgs <- chromosome_data %>%
  group_by(genotype, stage, chromosome) %>%
  summarise_all(mean)

for(variable in variable_name)
{

  p1 <- ggplot(chromosome_data, aes_(as.name(variable)), aes(group = genotype)) +
    geom_density(alpha = 0.8,aes(group = genotype)) +
    facet_wrap(~stage*chromosome) +
    theme_light() +
    geom_histogram(aes(y = ..density.., fill = genotype), alpha = 0.8,position = "dodge") 
  
  
  p2 <- ggplot(chromosome_data, aes_(y = as.name(variable), x = as.name('chromosome'), fill = as.name('genotype'))) +
      geom_boxplot() +
      theme_light() +
      facet_wrap(~stage) +
      stat_compare_means( label ="p.signif")
  
  p3 <- ggplot(chromosome_data, aes_(y = as.name(variable), x = as.name('genotype'), fill = as.name('stage'))) +
    geom_boxplot() +
    theme_light() +
    facet_wrap(~chromosome) +
    stat_compare_means( label ="p.signif") +
    scale_fill_brewer(palette = 'Set2')
  
  
  plot_list <- list(p1,p2,p3)
  
  ggsave(paste0('Plots/',variable, '.pdf'), marrangeGrob(grobs = plot_list, nrow = 1, ncol = 1))
  
}