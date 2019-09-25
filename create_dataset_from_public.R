unified_dataset_creator_only_x <- function(file_path, genotype, stage)
{
  chromosome_order <- c('chr_X')
  
  file_names <- list.files(file_path, full.names = T) 
  
  dendrites <- grep(file_names, pattern = 'Alignment_Dendrite', fixed = T, value = T)
  segments <- grep(file_names, pattern = 'Filament', fixed = T, value = T)
  
  nu_segments <- grep(segments, pattern = 'Dendrite_Branches|Sholl_Intersections', value = T, invert = T)
  
  segment_data <- nu_segments %>%
    map(read_csv,skip = 3) %>%
    lapply(function(x){x[,c(1, ncol(x)-1)]}) %>%
    reduce(full_join, by = c('ID' = 'ID')) 
  
  
  dendrite_data <- dendrites %>%
    map(read_csv,skip = 3) %>%
    lapply(function(x){x[,c(1, ncol(x)-2,ncol(x)-1)]}) %>%
    reduce(full_join, by = c('ID' = 'ID'))
  
  
  #get file name to then merge dendrites and filaments, it can come from any of the output files
  pic_name <- str_split(basename(dendrites[1]), pattern = '-')[[1]][1]
  
  all_segment <- segment_data %>%
    add_column(genotype = genotype,
               stage = stage, 
               chromosome = chromosome_order,
               file = pic_name) 
  
  all_dendrite <- dendrite_data %>%
    add_column(genotype = genotype,
               stage = stage,
               file = pic_name) %>%
    select(-contains('.x.')) %>%
    select(-contains('.y'))
  
  
  return(list(all_segment, all_dendrite))
  
}

library(tidyverse)

all_data_paths <- list.dirs(path = 'X:/Amalia/20190925_Imaris/', recursive = T) %>%
  as_data_frame() %>%
  filter(str_detect(value, pattern = 'Channel')) %>%
  mutate(stage = case_when(
    str_detect(value, 'EP') ~ 'Early Pachytene',
    str_detect(value, 'LP') ~ 'Late Pachytene'
  )) %>%
  mutate(base = basename(value)) %>%
  mutate(genotype = str_split(base, '-')) %>%
  unnest() %>%
  filter(!str_detect(genotype,'Channel')) %>%
  mutate(genotype = str_sub(genotype, start = 1, end = -5)) %>%
  mutate(genotype = case_when(
    str_detect(genotype, 'coh4') ~ str_sub(genotype, start = 1, end = -5),
    str_detect(genotype, 'cosa1_') ~ str_sub(genotype, start = 1, end = -6),
    TRUE ~ str_sub(genotype, start = 4)
  ))

data_list <- list()
for(i in 1:nrow(all_data_paths))
{
  data_list[[i]] <- unified_dataset_creator_only_x(file_path = all_data_paths$value[i],
                                 genotype = all_data_paths$genotype[i],
                                 stage = all_data_paths$stage[i])
}
  
filament <- lapply(data_list, function(x){x[[1]]}) %>%
  bind_rows() %>%
  write_csv('all_data_from_public.csv')

dendrite <- lapply(data_list, function(x){x[[2]]}) %>%
  bind_rows() %>%
  write_csv('all_dendrite_from_public.csv')
