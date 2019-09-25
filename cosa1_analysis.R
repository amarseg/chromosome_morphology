library(tidyverse)
library(ggpubr)

all_data <- read_csv('all_data_from_public.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

ggplot(filter(all_data, str_detect(genotype,'coh4')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means()
ggsave('coh4.pdf')

ggplot(filter(all_data, str_detect(genotype,'SMC1')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means()
ggsave('smc1.pdf')

ggplot(filter(all_data, str_detect(genotype,'rec8')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means()
ggsave('rec8.pdf')


ggplot(filter(all_data, str_detect(genotype,'cosa|spo')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means(ref.group = 'cosa1')
ggsave('cosa1_spo11.pdf')

summary_data <- all_data %>%
  group_by(genotype, stage) %>%
  summarise(n = n(), mean = mean(`Filament Length (sum)`), sd = sd(`Filament Length (sum)`)) %>%
  write_csv('summary_cosa1_exps.csv')

#####################Dendrites###################


chromosome_data <- read_csv('all_data_from_public.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))
dendrites_data <- read_csv('all_dendrite_from_public.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))


annot_dendrites <- inner_join(dendrites_data, chromosome_data, by= c('file' = 'file', 
                                                                     'FilamentID.x' = 'ID',
                                                                     'stage' = 'stage',
                                                                     'genotype' = 'genotype')) %>%
  group_by(file, FilamentID, stage, genotype) %>%
  mutate(segment_n = 1:n()) %>%
  filter(stage == 'Late Pachytene') %>%
  mutate(segment_ratio = `Dendrite Length`/`Filament Length (sum)`)


ggplot(filter(annot_dendrites, str_detect(genotype, 'SMC')), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.5)

ggplot(filter(annot_dendrites, str_detect(genotype, 'SMC')), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y=(..count..)/sum(..count..)),position = position_dodge()) 
ggsave('CO_SMC.pdf')


ggplot(filter(annot_dendrites, str_detect(genotype, 'coh')), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.5)

ggplot(filter(annot_dendrites, str_detect(genotype, 'coh')), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y=(..count..)/sum(..count..)),position = position_dodge()) 
ggsave('CO_coh.pdf')


ggplot(filter(annot_dendrites, str_detect(genotype, 'rec8')), aes(x = segment_ratio, fill = genotype))+
  geom_density(alpha = 0.5)

ggplot(filter(annot_dendrites, str_detect(genotype, 'rec8')), aes(x = segment_ratio, fill = genotype))+
  geom_histogram(aes(y=(..count..)/sum(..count..)),position = position_dodge()) 
ggsave('CO_rec.pdf')
