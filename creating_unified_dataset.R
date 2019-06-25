#######Function to read data and put it in a large table########
library('tidyverse')



unified_dataset_creator <- function(file_path, genotype, stage)
{
  chromosome_order <- c('chr_X','chr_III','chr_V')
  
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
  
  chromosome_names <- rep(chromosome_order,nrow(segment_data)/3)
  
  #get file name to then merge dendrites and filaments, it can come from any of the output files
  pic_name <- str_split(basename(dendrites[1]), pattern = '-')[[1]][1]
  
  all_segment <- segment_data %>%
    add_column(genotype = genotype,
               stage = stage, 
               chromosome = chromosome_names,
               file = pic_name) 
    
  all_dendrite <- dendrite_data %>%
    add_column(genotype = genotype,
               stage = stage,
               file = pic_name) %>%
    select(-contains('.x.')) %>%
    select(-contains('.y'))
    
  
  return(list(all_segment, all_dendrite))

}


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
wt_EP_1 <- unified_dataset_creator('DATASET_2/Wild Type/Early Pachytene/EP_MD_4_TIRF- Filtered_Channel Alignment_Statistics/',
                                genotype = 'wt',
                                stage = 'Early Pachytene')

wt_EP_2 <- unified_dataset_creator('DATASET_2/Wild Type/Early Pachytene/MP-LP_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wt',
                                   stage = 'Early Pachytene')

wapl_EP_1 <- unified_dataset_creator('DATASET_2/Wapl-1/Early Pachytene/wapl_EP_LP3_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wapl',
                                   stage = 'Early Pachytene')

wapl_EP_2 <- unified_dataset_creator('DATASET_2/Wapl-1/Early Pachytene/wapl_EP_LP5_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wapl',
                                   stage = 'Early Pachytene')

wt_LP_1 <- unified_dataset_creator('DATASET_2/Wild Type/Late Pachytene/LP_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wt',
                                   stage = 'Late Pachytene')

wt_LP_2 <- unified_dataset_creator('DATASET_2/Wild Type/Late Pachytene/LP2_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wt',
                                   stage = 'Late Pachytene')

wt_LP_3 <- unified_dataset_creator('DATASET_2/Wild Type/Late Pachytene/LP3_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wt',
                                   stage = 'Late Pachytene')

wt_LP_4 <- unified_dataset_creator('DATASET_2/Wild Type/Late Pachytene/LP9_TIRF- Filtered_Channel Alignment_Statistics/',
                                   genotype = 'wt',
                                   stage = 'Late Pachytene')


wapl_LP_1 <- unified_dataset_creator('DATASET_2/Wapl-1/Late Pachytene/wapl_EP_LP3_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'wapl',
                                     stage = 'Late Pachytene')

wapl_LP_2 <- unified_dataset_creator('DATASET_2/Wapl-1/Late Pachytene/wapl_EP_LP4_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'wapl',
                                     stage = 'Late Pachytene')
wapl_LP_3 <- unified_dataset_creator('DATASET_2/Wapl-1/Late Pachytene/wapl_EP_LP5_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'wapl',
                                     stage = 'Late Pachytene')
wapl_LP_4 <- unified_dataset_creator('DATASET_2/Wapl-1/Late Pachytene/wapl_EP2_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'wapl',
                                     stage = 'Late Pachytene')

wapl_LP_5 <- unified_dataset_creator('DATASET_2/Wapl-1/Late Pachytene/wapl_LP5_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'wapl',
                                     stage = 'Late Pachytene')

coh4DG_Auxin4_1_2 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/EP/coh4DG_Auxin4_1_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                     genotype = 'auxin_coh4',
                                     stage = 'Early Pachytene')


coh4DG_Auxin4_2_2 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/EP/coh4DG_Auxin4_2_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Early Pachytene')

coh4DG_Auxin4_3_2 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/EP/coh4DG_Auxin4_3_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Early Pachytene')

coh4DG_Auxin4_6_2 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/EP/coh4DG_Auxin4_6_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Early Pachytene')


coh4DG_Auxin4_1_1 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/LP/coh4DG_Auxin4_1_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Late Pachytene')


coh4DG_Auxin4_2_1 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/LP/coh4DG_Auxin4_2_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Late Pachytene')

coh4DG_Auxin4_3_1 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/LP/coh4DG_Auxin4_3_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Late Pachytene')

coh4DG_Auxin4_6_1 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/LP/coh4DG_Auxin4_6_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Late Pachytene')

coh4DG_Auxin4_8_1 <- unified_dataset_creator_only_x('coh4_DG/Auxin_coh4_DG/LP/coh4DG_Auxin4_8_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'auxin_coh4',
                                                    stage = 'Late Pachytene')


coh4DG_control_1_2 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/EP/coh4DG_control_1_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'control_coh4',
                                                    stage = 'Early Pachytene')


coh4DG_control_3_3 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/EP/coh4DG_control_3_3_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'control_coh4',
                                                    stage = 'Early Pachytene')

coh4DG_control_4_2 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/EP/coh4DG_control_4_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'control_coh4',
                                                    stage = 'Early Pachytene')

coh4DG_control_5_2 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/EP/coh4DG_control_5_2_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'control_coh4',
                                                    stage = 'Early Pachytene')

coh4DG_control_1_1 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/LP/coh4DG_control_1_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                    genotype = 'control_coh4',
                                                    stage = 'Late Pachytene')
coh4DG_control_3_1 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/LP/coh4DG_control_3_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                     genotype = 'control_coh4',
                                                     stage = 'Late Pachytene')
coh4DG_control_4_1 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/LP/coh4DG_control_4_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                     genotype = 'control_coh4',
                                                     stage = 'Late Pachytene')
coh4DG_control_5_1 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/LP/coh4DG_control_5_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                     genotype = 'control_coh4',
                                                     stage = 'Late Pachytene')

coh4DG_control_6_1 <- unified_dataset_creator_only_x('coh4_DG/control_coh4_DG/LP/coh4DG_control_6_1_TIRF- Filtered_Channel Alignment_Statistics/',
                                                     genotype = 'control_coh4',
                                                     stage = 'Late Pachytene')

filament_dataset <- bind_rows(wt_EP_1[[1]],wt_EP_2[[1]],
                              wapl_EP_1[[1]], wapl_EP_2[[1]], 
                              wt_LP_1[[1]],wt_LP_2[[1]],wt_LP_3[[1]],wt_LP_4[[1]], 
                              wapl_LP_1[[1]],wapl_LP_2[[1]],wapl_LP_3[[1]],wapl_LP_4[[1]],wapl_LP_5[[1]]) %>%
  write_csv('all_filaments_2.csv')


dendrite_dataset <- bind_rows(wt_EP_1[[2]],wt_EP_2[[2]],
                              wapl_EP_1[[2]], wapl_EP_2[[2]], 
                              wt_LP_1[[2]],wt_LP_2[[2]],wt_LP_3[[2]],wt_LP_4[[2]], 
                              wapl_LP_1[[2]],wapl_LP_2[[2]],wapl_LP_3[[2]],wapl_LP_4[[2]],wapl_LP_5[[2]]) %>%
  write_csv('all_dendrites_2.csv')


cos_filament_dataset <- bind_rows(coh4DG_Auxin4_1_2[[1]], coh4DG_Auxin4_2_2[[1]],coh4DG_Auxin4_3_2[[1]],coh4DG_Auxin4_6_2[[1]],
                                  coh4DG_Auxin4_1_1[[1]],coh4DG_Auxin4_2_1[[1]],coh4DG_Auxin4_3_1[[1]],coh4DG_Auxin4_6_1[[1]],
                                  coh4DG_Auxin4_8_1[[1]],
                                  coh4DG_control_1_2[[1]],coh4DG_control_3_3[[1]],coh4DG_control_4_2[[1]],
                                  coh4DG_control_5_2[[1]],coh4DG_control_1_1[[1]],coh4DG_control_3_1[[1]],
                                  coh4DG_control_4_1[[1]],coh4DG_control_5_1[[1]],coh4DG_control_6_1[[1]]) %>%
  write_csv('cos_filaments.csv')


cos_dendrite_dataset <- bind_rows(coh4DG_Auxin4_1_2[[2]], coh4DG_Auxin4_2_2[[2]],coh4DG_Auxin4_3_2[[2]],coh4DG_Auxin4_6_2[[2]],
                                  coh4DG_Auxin4_1_1[[2]],coh4DG_Auxin4_2_1[[2]],coh4DG_Auxin4_3_1[[2]],coh4DG_Auxin4_6_1[[2]],
                                  coh4DG_Auxin4_8_1[[2]],
                                  coh4DG_control_1_2[[2]],coh4DG_control_3_3[[2]],coh4DG_control_4_2[[2]],
                                  coh4DG_control_5_2[[2]],coh4DG_control_1_1[[2]],coh4DG_control_3_1[[2]],
                                  coh4DG_control_4_1[[2]],coh4DG_control_5_1[[2]],coh4DG_control_6_1[[2]]) %>%
write_csv('cos_dendrites.csv')





