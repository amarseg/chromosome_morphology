rm(list = ls())
library(tidyverse)
library(ggpubr)
library(gridExtra)

dendrites_data <- read_csv('all_dendrites.csv') 
chromosome_data <- read_csv('all_filaments.csv') 


annot_dendrites <- inner_join(dendrites_data, chromosome_data, by= c('file' = 'file', 
                                                                    'FilamentID' = 'ID',
                                                                    'stage' = 'stage',
                                                                    'genotype' = 'genotype')) %>%
  group_by(file, FilamentID, stage, genotype) %>%
  mutate(segment_n = 1:n()) %>%
  filter(stage == 'Late Pachytene') %>%
  mutate(segment_ratio = `Dendrite Length`/`Filament Length (sum)`)

ggplot(annot_dendrites, aes(y = segment_ratio, x = chromosome, fill = genotype))+
  geom_boxplot() +
  facet_wrap(~segment_n) +
  stat_compare_means( label ="p.signif")


ggplot(filter(annot_dendrites, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y = ..density..),position = "dodge") +
  facet_wrap(~chromosome, nrow = 3, ncol = 1) 

ggplot(filter(annot_dendrites, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.75) +
  facet_wrap(~chromosome, nrow = 3, ncol = 1) 
