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

binned_things <- annot_dendrites %>%
  group_by(genotype, chromosome) %>%
  mutate(bin = ntile(`Filament Length (sum)`,2)) %>%
  ungroup()

ggplot(binned_things, aes(y = `Filament Length (sum)`, x = chromosome, fill = genotype)) +
  geom_boxplot() +
  theme_light() +
  facet_wrap(~bin) 

ggplot(binned_things, aes(x = segment_ratio, fill = genotype)) +
  geom_histogram(aes(y = ..density..),position = "dodge") +
  facet_wrap(~bin*chromosome) 

ggplot(binned_things, aes(x = segment_ratio, fill = as.factor(bin))) +
  geom_density(aes(y = ..density..),position = "dodge", alpha = 0.75) +
  facet_wrap(~genotype*chromosome) 
