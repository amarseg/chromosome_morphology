library(tidyverse)
library(ggpubr)

all_data <- read_csv('all_data_from_public.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

ggplot(filter(all_data, str_detect(genotype,'coh4')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  stat_compare_means()
ggsave('coh4.pdf')

ggplot(filter(all_data, str_detect(genotype,'SMC1')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  stat_compare_means()
ggsave('smc1.pdf')

ggplot(filter(all_data, str_detect(genotype,'rec8')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  stat_compare_means()
ggsave('rec8.pdf')


ggplot(filter(all_data, str_detect(genotype,'cosa|spo')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  stat_compare_means(ref.group = 'cosa1')
ggsave('cosa1_spo11.pdf')

summary_data <- all_data %>%
  group_by(genotype, stage) %>%
  summarise(n = n(), mean = mean(`Filament Length (sum)`), sd = sd(`Filament Length (sum)`)) %>%
  write_csv('summary_cosa1_exps.csv')
