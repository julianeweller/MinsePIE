# 2021-11-08
# This script takes the read count tables generated in read_match_pegRNAs script, merges replicates and adds additional filtering and normalization steps resulting in Table_S2_insertion_frequencies.csv
library(tidyverse)

# ==== Loading files ====
# Run this function with a table with all file names as input, a path to the directory with raw read count tables from the Read_match_pegRNA script, and a string indicating the column name of the sample column
# Example: Gall_prc <- read_files(prc_Gall_file_names, "/path/to/directory/", "sample")
read_files <- function(files, path, indicator) {
  files <- mutate(files, path = paste0(path, file))
  data = list()
  
  for(i in 1:nrow(files)) {
    data[[i]] = read_tsv(files$path[i])
  }
  
  names <- vector("character", length = length(data))
  for(i in 1:length(data)) {
    names[i] = (data[[i]][[indicator]][1])
  }
  
  names(data) = names
  
  return(data)
}

# ==== Functions ====
# Function to merge replicates and generate ratios from pegRNA and target read counts
# Run function with data frames from the table generated from the read_files function, providing three tables for the three pegRNA replicates and three tables for the three target replicates
# Example: FANCF_293T <- combine_replicates(Gall_prc$`45-1-1-pegRNA`,Gall_prc$`45-1-2-pegRNA`, Gall_prc$`45-1-3-pegRNA`, `45-2-1-target`, `45-2-2-target`, `45-2-3-target`) %>% mutate(experiment = "FANCF_293T_PE2_3")
combine_replicates <- function(R1_peg, R2_peg, R3_peg, R1_target, R2_target, R3_target) {
  # Merging pegRNAs
  pegs <- full_join(R1_peg, R2_peg, by = "sequence_original") %>%
    full_join(R3_peg, by = "sequence_original") %>%
    select(name = name.x, sequence_original, count_peg_R1 = count_norm.x, count_peg_R2 = count_norm.y, count_peg_R3 = count_norm)
  # Merging targets
  targets <- full_join(R1_target, R2_target, by = "sequence_original") %>%
    full_join(R3_target, by = "sequence_original") %>%
    select(sequence_original, count_target_R1 = count_norm.x, count_target_R2 = count_norm.y, count_target_R3 = count_norm)
  # Calculating insertion frequencies and augmenting the resulting data frame
  ratio <- full_join(targets, pegs, by = "sequence_original") %>% 
    mutate(percIns_R1 = count_target_R1/count_peg_R1*100, percIns_R2 = count_target_R2/count_peg_R2*100, percIns_R3 = count_target_R3/count_peg_R3*100, 
           percIns = (percIns_R1 + percIns_R2 + percIns_R3)/3)
}

# Combining data in one long data frame for filtering
combined_all <- bind_rows(FANCF_293T, HEK3_293T, FANCF_HAP1, HEK3_HAP1, HEK3_293T_FeLV, CLYBL_293T, EMX1_293T, EMX1_nicking_293T) %>% 
  mutate(length = nchar(sequence_original)) %>%
  left_join(distinct(library_VF), by = "sequence_original") %>% # Provide a file with Vienna fold values for each pegRNA
  left_join(distinct(library_annotation), by = "sequence_original")  %>% # Provide an optional file with additional library annotations
  left_join(tibble(length = seq(0:69), bin = c(rep(seq(5, 30, by = 5), each = 5), rep(seq(40, 70, by = 10), each = 10))), by = "length") %>% # Generating length bins
  # This step removes any duplicates from reverse complements
  group_by(sequence_original, experiment) %>%
  filter(row_number() == 1)

# Filtering and normalization
combined_all <-  combined_all %>%
  # First require at least 10 reads for each pegRNA in each replicate
  filter(count_peg_R1 > 10, count_peg_R2 > 10, count_peg_R3 > 10) %>%
  # Generate a Vienna fold column that takes the matching VF fold value for each insert
  mutate(VF_full = ifelse(str_detect(experiment, "CLYBL"), VF_full_CLYBL, ifelse(
    str_detect(experiment, "EMX1"), VF_full_EMX1, ifelse(
      str_detect(experiment, "FANCF"), VF_full_FANCF, VF_full_HEK3))),
    # Generating parameters
    percGC = (str_count(sequence_original, "G") + str_count(sequence_original, "C"))/nchar(sequence_original) *100,
    percA = str_count(sequence_original, "A")/nchar(sequence_original) *100,
    percC = str_count(sequence_original, "C")/nchar(sequence_original) *100,
    percT = str_count(sequence_original, "T")/nchar(sequence_original) *100,
    percG = str_count(sequence_original, "G")/nchar(sequence_original) *100,
    Trun = str_detect(sequence_original, "TTTT"), Arun = str_detect(sequence_original, "AAAA"), Crun = str_detect(sequence_original, "CCCC"), Grun = str_detect(sequence_original, "GGGG")) %>%
  # Data normalization
  group_by(experiment) %>%
  mutate(percIns_z = (percIns-mean(percIns, na.rm = T))/sd(percIns, na.rm = T), percIns_norm = percIns/median(percIns, na.rm = T)) %>%
  group_by(experiment, bin) %>%
  mutate(len_res = percIns / median(percIns, na.rm = T), VF_full_norm = VF_full/median(VF_full, na.rm = T), VF_insert_norm = VF_insert/(median(VF_insert, na.rm = T)+0.0001)) %>%
  ungroup()

# The combined_all tibble generated in this script corresponds to Table_S2_insertion_frequencies.csv in the supplementary information