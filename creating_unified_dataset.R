#######Function to read data and put it in a large table########
library('tidyverse')
chromosome_order <- c('chr_X','chr_III','chr_V')



#####Test case################
file_path <- 'Amalia/Wild Type/Early Pachytene/EP_MD_4_TIRF- Filtered_Channel Alignment_Statistics/'
file_names <- list.files(file_path, full.names = T)

chromosomes <- grep(file_names, pattern = 'Alignment_Dendrite', fixed = T, value = T)
segments <- grep(file_names, pattern = 'Filament', fixed = T, value = T) %>%
  select(-(str_detect()))
parameter_name <- str_sub(basename(segments), start = 42, end = nchar(basename(segments)) - 4)

chromosome_data <- segments %>%
  map(read_csv, skip = 2) %>%
  
  reduce(cbind)
