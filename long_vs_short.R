library(tidyverse)
library(ggpubr)
library(gridExtra)
library(wesanderson)

lp_data <- read_csv('processed_LP.csv') %>%
  select( Genotype, Chromosome, Avg_Long_angle_1, Avg_Short_angle_1,Avg_Around_angle_1,
          Avg_Long_angle_2, Avg_Short_angle_2,Avg_Around_angle_2) %>%
  gather(key = 'key', value = 'measurement', -Chromosome, -Genotype) %>%
  filter(str_detect(key, 'angle_1'))
  

ggplot(lp_data, aes(x = key, fill = Genotype, y = measurement)) +
  geom_boxplot() + 
  facet_wrap(~Chromosome) +
  stat_compare_means(label ="p.signif") +
  scale_fill_brewer(palette = 'Accent') 
ggsave('fig_3d/short_vs_long.png')

ggplot(lp_data, aes(x = Chromosome, fill = key, y = measurement)) +
  geom_boxplot() +
  facet_wrap(~Genotype) +
  stat_compare_means( label ="p.signif", method = 'anova') +
  scale_fill_brewer(palette = 'Accent') +
  title('Anova')
ggsave('fig_3d/short_vs_long_anova.png')

ggplot(lp_data, aes(x = Chromosome, fill = Genotype, y = measurement)) +
  geom_boxplot() +
  facet_wrap(~key) +
  stat_compare_means(label ="p.signif", method = 'anova') +
  title('Anova')
ggsave('fig_3d/short_vs_long_genotype_anova.png')

ggplot(lp_data, aes(x = Chromosome, fill = Genotype, y = measurement)) +
  geom_boxplot() +
  facet_wrap(~key) +
  stat_compare_means(label ="p.signif")
ggsave('fig_3d/short_vs_long_genotype.png')