library(tidyverse)
library(ggpubr)
library(gridExtra)
library(lawstat)

spo_fil <- read_csv('Spo_filaments.csv')
rest_data <- read_csv('all_filaments_2.csv') %>%
  filter(chromosome == 'chr_X')

all_data <- bind_rows(spo_fil, rest_data)

ggplot(all_data, aes(y = `Filament Length (sum)`, x = genotype)) +
  geom_boxplot() +
  facet_wrap(~stage)+
  stat_compare_means( label ="p.signif", ref.group = 'wt')

ggplot(all_data, aes(y = `Filament Length (sum)`, x = stage)) +
  geom_boxplot() +
  facet_wrap(~genotype)

only_wt <- read_csv('all_filaments_2.csv') %>%
  filter(chromosome == 'chr_X' & genotype == 'wt')

ggplot(filter(all_data, genotype != 'wapl'), aes(x = `Filament Length (sum)`, fill = genotype)) +
  geom_histogram(position = position_dodge()) +
  facet_wrap(~stage)
