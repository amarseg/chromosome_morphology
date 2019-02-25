rm(list = ls())
library(tidyverse)
library(ggpubr)
library(gridExtra)

dendrites_data <- read_csv('all_dendrites_2.csv') 
chromosome_data <- read_csv('all_filaments_2.csv') 


annot_dendrites <- inner_join(dendrites_data, chromosome_data, by= c('file' = 'file', 
                                                                    'FilamentID.x' = 'ID',
                                                                    'stage' = 'stage',
                                                                    'genotype' = 'genotype')) %>%
  group_by(file, FilamentID.x, stage, genotype) %>%
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

ggplot(filter(annot_dendrites, segment_n ==1), aes(x = segment_ratio, y = `Filament Length (sum)`,
                                                   color = chromosome)) +
  geom_point()

test_data <- chromosome_data %>%
  mutate(edges_per_length = `Filament No. Edges`/`Filament Length (sum)`)

ggplot(test_data, aes(y = edges_per_length, x = genotype, fill = stage)) +
  geom_boxplot() +
  theme_light() +
  facet_wrap(~chromosome) +
  stat_compare_means( label ="p.signif") +
  scale_fill_brewer(palette = 'Set2')

ggplot(test_data, aes(y = edges_per_length, x = chromosome, fill = genotype)) +
  geom_boxplot() +
  theme_light() +
  facet_wrap(~stage) +
  stat_compare_means( label ="p.signif") 
