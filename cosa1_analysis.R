library(tidyverse)
library(ggpubr)

all_data <- read_csv('all_data_from_public.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

all_data2 <- read_csv('all_data_from_public2.csv') %>%
  mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

all_data <- bind_rows(all_data, all_data2)

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
  stat_compare_means(ref.group = '0Auxin5_rec8GFPDG_',label = "p.signif")
ggsave('rec8.pdf')


ggplot(filter(all_data, str_detect(genotype,'cosa|spo')),
       aes( x= genotype, y = `Filament Length (sum)`)) +
  geom_boxplot() +
  facet_wrap(~stage) +
  stat_compare_means(ref.group = 'cosa1', label = "p.signif")
ggsave('cosa1_spo11.pdf')

summary_data <- all_data %>%
  group_by(genotype, stage) %>%
  summarise(n = n(), mean = mean(`Filament Length (sum)`), sd = sd(`Filament Length (sum)`)) %>%
  write_csv('summary_cosa1_exps.csv')

#####################Dendrites###################


# chromosome_data <- read_csv('all_data_from_public.csv') %>%
#   mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))
# dendrites_data <- read_csv('all_dendrite_from_public.csv') %>%
#   mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

chromosome_data2 <- read_csv('all_data_from_public2.csv') %>%
   mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))
dendrites_data2 <- read_csv('all_dendrite_from_public2.csv') %>%
   mutate(genotype = str_replace(genotype, pattern = '_', replacement = ''))

chromosome_data <- bind_rows(chromosome_data, chromosome_data2)
dendrites_data <- bind_rows(dendrites_data, dendrites_data2)

annot_dendrites <- inner_join(dendrites_data2, chromosome_data2, by= c('file' = 'file', 
                                                                     'FilamentID.x' = 'ID',
                                                                     'stage' = 'stage',
                                                                     'genotype' = 'genotype')) %>%
  group_by(file, FilamentID.x, stage, genotype) %>%
  mutate(segment_n = 1:n()) %>%
  
  filter(stage == 'Late Pachytene') %>%
  mutate(ratio_fragment = `Dendrite Length`/`Filament Length (sum)`) %>%
  group_by(file, FilamentID.x, stage, genotype, chromosome) %>%
  summarise(min_fragment = min(ratio_fragment), mean_length = mean(`Filament Length (sum)`),
            diff_ratio = max(ratio_fragment) - min(ratio_fragment)) 

ggplot(filter(annot_dendrites, str_detect(genotype, 'SMC')), aes(x = min_fragment, fill = genotype))+
  geom_density(alpha = 0.5)

ggplot(filter(annot_dendrites, str_detect(genotype, 'SMC')), aes(x = min_ratio, fill = genotype))+
  geom_histogram(aes(y=(..count..)/sum(..count..)),position = position_dodge()) 
ggsave('CO_SMC.pdf')


ggplot(filter(annot_dendrites, str_detect(genotype, 'coh')), aes(x = min_ratio, fill = genotype))+
  geom_density(alpha = 0.5)

ggplot(filter(annot_dendrites, str_detect(genotype, 'coh')), aes(x = min_ratio, fill = genotype))+
  geom_histogram(aes(y=(..count..)/sum(..count..)),position = position_dodge()) 
ggsave('CO_coh.pdf')


ggplot(filter(annot_dendrites, str_detect(genotype, 'rec8')), aes(x = min_ratio, fill = genotype))+
  geom_density(alpha = 0.5)

library(scales)
ggplot(filter(annot_dendrites, str_detect(genotype, 'rec8')), aes(x = min_ratio, fill = genotype, group = genotype))+
  geom_histogram(position = position_dodge()) +
  scale_fill_brewer(palette = 'Dark2') +
  facet_wrap(~genotype, scales = 'free_y')

ggsave('CO_rec.pdf')
