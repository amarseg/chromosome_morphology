library(tidyverse)
library(ggpubr)
library(gridExtra)

###################
##Plot Length data

chromosome_data <- read_csv('all_filaments_2.csv') %>%
  select(genotype, stage, chromosome, `Filament Length (sum)`, ID, file) %>%
  filter(chromosome == 'chr_X')

degron_data <- read_csv('cos_filaments.csv') %>%
  select(genotype, stage, chromosome, `Filament Length (sum)`, ID, file)

ggplot(degron_data, aes(x = genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means(method = )

ggplot(degron_data, aes(x = stage, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~genotype) +
  stat_compare_means()

all_data <- bind_rows(chromosome_data, degron_data)

ggplot(all_data, aes(x = genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means(ref.group = 'wt', label = "p.format") 
ggsave('Plots/degron/length_vs_genotype.pdf', width = 10)

ggplot(all_data, aes(x = stage, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~genotype) +
  stat_compare_means(ref.group = 'wt')
ggsave('Plots/degron/length_vs_stage.pdf')


ggplot(filter(all_data, stage == 'Early Pachytene'), aes(`Filament Length (sum)`, fill = genotype)) +
  geom_density(alpha = 0.75) 
ggsave('Plots/degron/overlap_length_distribution_EP.pdf')

ggplot(filter(all_data, stage == 'Late Pachytene'), aes(`Filament Length (sum)`, fill = genotype)) +
  geom_density(alpha = 0.75) 
ggsave('Plots/degron/overlap_length_distribution_LP.pdf')
#####################
#Plot crossover data

dendrites_data <- read_csv('all_dendrites_2.csv') 


annot_dendrites <- inner_join(dendrites_data, chromosome_data, by= c('file' = 'file', 
                                                                     'FilamentID.x' = 'ID',
                                                                     'stage' = 'stage',
                                                                     'genotype' = 'genotype')) %>%
  group_by(file, FilamentID.x, stage, genotype) %>%
  mutate(segment_n = 1:n()) %>%
  filter(stage == 'Late Pachytene') %>%
  mutate(segment_ratio = `Dendrite Length`/`Filament Length (sum)`)

degron_dendrites <- read_csv('cos_dendrites.csv')

annot_degron <- inner_join(degron_dendrites, degron_data, by= c('file' = 'file', 
                                                                  'FilamentID.x' = 'ID',
                                                                  'stage' = 'stage',
                                                                  'genotype' = 'genotype')) %>%
  group_by(file, FilamentID.x, stage, genotype) %>%
  mutate(segment_n = 1:n()) %>%
  filter(stage == 'Late Pachytene') %>%
  mutate(segment_ratio = `Dendrite Length`/`Filament Length (sum)`)


ggplot(filter(annot_degron, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y = ..density..),position = "dodge") +
  facet_wrap(~chromosome, nrow = 3, ncol = 1) 

ggplot(filter(annot_degron, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.75) 


annot_all <- bind_rows(annot_dendrites, annot_degron)

ggplot(filter(annot_all, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y = ..density..),position = "dodge") +
  facet_wrap(~genotype) 
ggsave('Plots/degron/histogram_crossover.pdf')

ggplot(filter(annot_all, segment_n == 1), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.75) +
  facet_wrap(~genotype)
ggsave('Plots/degron/density_crossover.pdf')

summary_stats <- all_data %>%
  group_by(stage, genotype) %>%
  summarise(avg_length = mean(`Filament Length (sum)`), sd_length = sd(`Filament Length (sum)`)) %>%
  write_csv('summary_chromosome_stats.csv')
