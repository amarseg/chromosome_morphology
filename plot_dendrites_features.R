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
  facet_wrap(~genotype*chromosome, nrow = 2, ncol = 3) 

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

chr_x <- chromosome_data %>%
  filter(chromosome == 'chr_X') %>%
  arrange(ID) %>%
  select(ID, file, stage, `Filament Length (sum)`, genotype)

chr_iii <- chromosome_data %>%
  filter(chromosome == 'chr_III')   %>%
  arrange(ID) %>%
  select(ID, file, stage, `Filament Length (sum)`)

chr_v <- chromosome_data %>%
  filter(chromosome == 'chr_V')   %>%
  arrange(ID) %>%
  select(ID, file, stage, `Filament Length (sum)`)

t <- bind_cols(chr_x, chr_iii, chr_v)


p1<-ggplot(t, aes(x = `Filament Length (sum)`, y =  `Filament Length (sum)1`, group = genotype)) +
  geom_point(aes(color = stage)) +
  geom_smooth(method = 'lm', colour = 'black') +
  xlab('Chromosome X length') +
  ylab('Chromosome III length') +
  scale_colour_brewer(palette = 'Set2') +
  theme_bw() +
  facet_wrap(~genotype) +
  stat_cor(output.type = 'text')

p2 <-ggplot(t, aes(x = `Filament Length (sum)1`, y =  `Filament Length (sum)2`)) +
  geom_point(aes(color = stage)) +
  geom_smooth(method = 'lm', colour = 'black') +
  xlab('Chromosome III length') +
  ylab('Chromosome V length') +
  scale_colour_brewer(palette = 'Set2') +
  theme_bw() +
  facet_wrap(~genotype)  +
  stat_cor(output.type = 'text')

p3 <- ggplot(t, aes(x = `Filament Length (sum)`, y =  `Filament Length (sum)2`)) +
  geom_point(aes(color = stage)) +
  geom_smooth(method = 'lm', colour = 'black') +
  xlab('Chromosome  X length') +
  ylab('Chromosome V length') +
  scale_colour_brewer(palette = 'Set2') +
  theme_bw() +
  facet_wrap(~genotype) +
  stat_cor(output.type = 'text')

grid.arrange(p1,p2,p3, ncol = 1, nrow = 3)


####################
#Select smalles fraction for plots and do t-test 

min_dendrites <- annot_dendrites %>%
  group_by(file, FilamentID.x, stage, genotype, chromosome) %>%
  summarise(min_fragment = min(segment_ratio), mean_length = mean(`Filament Length (sum)`))


ggplot(min_dendrites, aes(y = min_fragment, x = chromosome, fill = genotype))+
  geom_boxplot() +
  stat_compare_means( label ="p.signif")

ggsave('Plots/crossover_min_wt.pdf')

ggplot(min_dendrites, aes(x = min_fragment, fill = genotype))+
  geom_histogram(aes(y = ..density..),position = "dodge") +
  facet_wrap(~chromosome, nrow = 3, ncol = 1) 

ggplot(min_dendrites, aes(x = min_fragment, fill = genotype))+
  geom_density(alpha = 0.75) +
  facet_wrap(~chromosome, nrow = 2, ncol = 3) 

ggplot(min_dendrites, aes(x = min_fragment, y = mean_length)) +
  geom_point(aes(colour = genotype)) +
  facet_wrap(~chromosome) +
  geom_smooth(method = 'lm')

ggsave('Plots/relationship_length_CO.pdf')
