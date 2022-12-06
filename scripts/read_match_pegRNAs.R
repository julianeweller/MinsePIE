# 2022-11-11
# This script takes fastq sequencing files as input and generates read count tables with library insertions

# Required libraries
library(ShortRead)
library(stringr)
library(tidyverse)

# A function to generate reverse complement nucleotides
rc <- function(x) {toupper(spgs::reverseComplement(x))}

# specifying paths
workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies"
workdir = "/Users/jw38/OneDrive/Prime_editing_efficiencies"

# Loading in modified libraries to detect prime editing events or pegRNAs
library_HEK3 <- read_tsv(paste0(workdir, "/files/library_HEK3.tsv"))
library_HEK3_target <- read_tsv(paste0(workdir, "/files/library_HEK3_target.tsv"))
library_FANCF <- read_tsv(paste0(workdir, "/files/library_FANCF.tsv"))
library_FANCF_target <- read_tsv(paste0(workdir, "/files/library_FANCF_target.tsv"))
library_CLYBL <- read_tsv(paste0(workdir, "/files/library_CLYBL.tsv"))
library_CLYBL_target <- read_tsv(paste0(workdir, "/files/library_CLYBL_target.tsv"))
library_EMX1 <- read_tsv(paste0(workdir, "/files/library_EMX1.tsv"))
library_EMX1_target <- read_tsv(paste0(workdir, "/files/library_EMX1_target.tsv"))
library_HEK3_RT <- mutate(library_HEK3, sequence = paste0("TTTCCTCTGCCATCA", rc(sequence_original), "CGTGCTCAGTCT"))
library_canadian_SP1 <- read_tsv(paste0(workdir, "/files/library_canadian_SP1.tsv"))
library_canadian_SP1_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP1_target.tsv"))
library_canadian_SP2 <- read_tsv(paste0(workdir, "/files/library_canadian_SP2.tsv"))
library_canadian_SP2_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP2_target.tsv"))
library_canadian_SP3 <- read_tsv(paste0(workdir, "/files/library_canadian_SP3.tsv"))
library_canadian_SP3_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP3_target.tsv"))
library_canadian_SP4 <- read_tsv(paste0(workdir, "/files/library_canadian_SP4.tsv"))
library_canadian_SP4_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP4_target.tsv"))
library_canadian_SP5 <- read_tsv(paste0(workdir, "/files/library_canadian_SP5.tsv"))
library_canadian_SP5_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP5_target.tsv"))
library_canadian_SP6 <- read_tsv(paste0(workdir, "/files/library_canadian_SP6.tsv"))
library_canadian_SP6_target <- read_tsv(paste0(workdir, "/files/library_canadian_SP6_target.tsv"))
library_barnacle <- read_tsv(paste0(workdir, "/files/library_barnacle.tsv"))
library_epeg <- read_tsv(paste0(workdir, "/files/library_epeg.tsv"))
library_epeg_target <- read_tsv(paste0(workdir, "/files/library_epeg_target.tsv"))




# A function to read and trim sequences and then detect pegRNAs / editing events
read_match_pegRNAs <- function(workdir, datapath, experiment, leftcons, rightcons, library) {
  
  # Paths
  path = paste0(workdir, datapath)
  
  # Extract the useful parts of the file name
  filename = sub("_S.*$", "", basename(path))
  print(filename)
  
  print("Reading fastq files")
  reads <- readFastq(path)
  reads <- as.character(reads@sread)
  # Trim reads
  reads_trim <- str_subset(reads, paste0(leftcons, "(.*?)", rightcons)) %>%
    str_extract(paste0(leftcons, "(.*?)", rightcons))
  
  perc_reads <- length(reads_trim)/length(reads)*100
  print(paste0("Percent of reads matched: ", perc_reads))

  print("Counting and matching pegRNAs")
  count_table <- tibble(sequence = reads_trim) %>%
    # First collapse reads
    group_by(sequence) %>%
    summarize(count = n()) %>%
    # Matching with the library
    full_join(library, by = "sequence") %>%
    # Identifying pegRNAs that were not detected in sequencing data
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    dplyr::select(name, sequence, sequence_original, count)
  
  perc_all <- sum(filter(count_table, !is.na(name))$count)/sum(count_table$count)*100
  perc_wt <- sum(filter(count_table, name == "WT")$count)/sum(count_table$count)*100
  perc_library <- sum(filter(count_table, !is.na(name) & name != "WT")$count)/sum(count_table$count)*100
  
  print(paste0("Percent of reads matched: ", perc_all))
  print(paste0("Percent of inserts matched: ", perc_library))
  
  # Some clean up. Fixing duplicatred sequences in the library. Normalization
  count_table <- filter(count_table, !is.na(name)) %>%
    group_by(sequence_original) %>%
    arrange(name) %>%
    summarise(name = name[1], sequence = sequence[1], count = count[1]/n()) %>%
    ungroup() %>%
    mutate(count_norm = count/sum(count)*1000000, sample = filename) %>%
    arrange(desc(count))
  
  perc_dropouts <- nrow(filter(count_table, count == 0))/nrow(count_table)*100
  n_dropouts <- nrow(filter(count_table, count == 0))
  avg_reads <- sum(filter(count_table, name != "WT")$count)/nrow(filter(count_table, name != "WT"))
  
  write_tsv(count_table, paste0(workdir, "/prc_data/", experiment, "_", filename, "_counts.txt"))
  write_tsv(tibble(perc_reads, perc_all, perc_library, perc_wt, perc_dropouts, n_dropouts, avg_reads, sum_counts = sum(count_table$count), sample = filename), paste0(workdir, "/prc_data/", experiment, "_", filename, "_metrics.txt"))
  
}

# works when there are several different constant regions
read_match_pegRNAs2 <- function(workdir, datapath, experiment, library) {
  
  # Paths
  path = paste0(workdir, datapath)
  
  # Extract the useful parts of the file name
  filename = sub("_S.*$", "", basename(path))
  print(filename)
  
  print("Reading fastq files")
  reads <- readFastq(path)
  reads <- as.character(reads@sread)
  
  leftcons = unique(library$leftcons)
  rightcons = unique(library$rightcons)
  print(leftcons)
  print(rightcons)
  
  # Trim reads
  reads_trim = vector(mode = "character")
  for(i in 1:length(leftcons)) {
    reads_trim_temp <- str_subset(reads, paste0(leftcons[i], "(.*?)", rightcons[i])) %>%
      str_extract(paste0(leftcons[i], "(.*?)", rightcons[i]))
    reads_trim = c(reads_trim, reads_trim_temp)
  }
  
  perc_reads <- length(reads_trim)/length(reads)*100
  print(paste0("Percent of reads matched: ", perc_reads))
  
  print("Counting and matching pegRNAs")
  count_table <- tibble(sequence = reads_trim) %>%
    # First collapse reads
    group_by(sequence) %>%
    summarize(count = n()) %>%
    # Matching with the library
    full_join(library, by = "sequence") %>%
    # Identifying pegRNAs that were not detected in sequencing data
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    dplyr::select(name, sequence, sequence_original, count)
  
  perc_all <- sum(filter(count_table, !is.na(name))$count)/sum(count_table$count)*100
  perc_wt <- sum(filter(count_table, name == "WT")$count)/sum(count_table$count)*100
  perc_library <- sum(filter(count_table, !is.na(name) & name != "WT")$count)/sum(count_table$count)*100
  
  print(paste0("Percent of reads matched: ", perc_all))
  print(paste0("Percent of inserts matched: ", perc_library))
  
  # Some clean up. Fixing duplicatred sequences in the library. Normalization
  count_table <- filter(count_table, !is.na(name)) %>%
    group_by(sequence) %>%
    arrange(name) %>%
    summarise(name = name[1], sequence_original = sequence_original[1], count = count[1]/n()) %>%
    ungroup() %>%
    mutate(count_norm = count/sum(count)*1000000, sample = filename) %>%
    arrange(desc(count))
  
  perc_dropouts <- nrow(filter(count_table, count == 0))/nrow(count_table)*100
  n_dropouts <- nrow(filter(count_table, count == 0))
  avg_reads <- sum(filter(count_table, name != "WT")$count)/nrow(filter(count_table, name != "WT"))
  
  write_tsv(count_table, paste0(workdir, "/prc_data/", experiment, "_", filename, "_counts.txt"))
  write_tsv(tibble(perc_reads, perc_all, perc_library, perc_wt, perc_dropouts, n_dropouts, avg_reads, sum_counts = sum(count_table$count), sample = filename), paste0(workdir, "/prc_data/", experiment, "_", filename, "_metrics.txt"))
  
}

# Processing of screen data
#45-1-4
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-4-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-4-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)

#45-1-5
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-5-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-5-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)

#45-1-6
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-6-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-6-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)


#45-1-13
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-13-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-13-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)
#45-1-14
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-14-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-14-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)
#45-1-15
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-15-target_S0_R1.fastq.gz",
                   experiment = "Goose", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210609_hiseq_walk-up_150_Goose-1/45-1-15-pegRNA_S0_R1.fastq.gz",
                   experiment = "Goose", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)



# HAP1

#45-2-7
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/FASTQ_Generation_2021-07-09_07_51_49Z-437178743/45-2-7-target_L001-ds.4736c22d65e949a9b9e2c59c603b188c/45-2-7-target_S7_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-7-pegRNA_S1_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)

#45-2-8
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-8-target_S8_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-8-pegRNA_S2_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)

#45-2-9
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-9-target_S9_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-9-pegRNA_S3_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)

#45-2-10
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-9-target_S9_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/MiSeq_walk_up_373-278374103/45-2-9-pegRNA_S3_L001_R1_001.fastq.gz",
                   experiment = "HAP1", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)


# FANCF 293T
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-2-1-target_S19_LRmerged.fastq.gz",
                   experiment = "Goose3", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-2-2-target_S20_LRmerged.fastq.gz",
                   experiment = "Goose3", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-2-3-target_S21_LRmerged_split1.fastq.gz",
                   experiment = "Goose3", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)

# EMX1 and CLYBL loci
#45-3-1
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-1-target_S10_L002_R2_001.fastq.gz",
                   experiment = "Goose3", "CAGACTGTCA", "CTAAGGACCT", library_CLYBL_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-1-pegRNA_S1_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "GACACAGGTCCTTAG", "TGACAGTCTGCACTT", library_CLYBL)

#45-3-2
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-2-target_S11_L001_R2_001.fastq.gz",
                   experiment = "Goose3", "CAGACTGTCA", "CTAAGGACCT", library_CLYBL_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-2-pegRNA_S2_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "GACACAGGTCCTTAG", "TGACAGTCTGCACTT", library_CLYBL)

#45-3-3
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-3-target_S12_L001_R2_001.fastq.gz",
                   experiment = "Goose3", "CAGACTGTCA", "CTAAGGACCT", library_CLYBL_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-3-pegRNA_S3_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "GACACAGGTCCTTAG", "TGACAGTCTGCACTT", library_CLYBL)


#45-3-4
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-4-target_S13_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-4-pegRNA_S4_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)

#45-3-5
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-5-target_S14_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-5-pegRNA_S5_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)

#45-3-6
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-6-target_S15_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-6-pegRNA_S6_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)



# With nicking guides
#45-3-4N
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-4N-target_S16_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-4N-pegRNA_S7_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)

#45-3-5N
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-5N-target_S17_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-5N-pegRNA_S8_L002_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)

#45-3-6N
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-6N-target_S18_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "AGCAGAAGAA", "GAAGGGCTCC", library_EMX1_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20210919_hiseq_walk-up_162_Goose-3/45-3-6N-pegRNA_S9_L001_R1_001.fastq.gz",
                   experiment = "Goose3", "TGATGGGAGCCCTTC", "TTCTTCTGCTCGGTT", library_EMX1)



#45-4-1
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-1-target_S16_RL-merged.fastq.gz",
                   experiment = "G4", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-1-pegRNA_S1_RL-merged.fastq.gz",
                   experiment = "G4", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)
#45-4-2
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-2-target_S17_RL-merged.fastq.gz",
                   experiment = "G4", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-2-pegRNA_S2_RL-merged.fastq.gz",
                   experiment = "G4", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)
#45-4-3
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-3-target_S18_RL-merged.fastq.gz",
                   experiment = "G4", "CTTCTGCAGC", "ACCTGGATCG", library_FANCF_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-3-pegRNA_S3_RL-merged.fastq.gz",
                   experiment = "G4", "AAAAGCGATCCAGGT", "GCTGCAGAAGGGATT", library_FANCF)

#45-4-4
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-4-target_S19_RL-merged.fastq.gz",
                   experiment = "G4", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-4-pegRNA_S4_RL-merged.fastq.gz",
                   experiment = "G4", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)
#45-4-5
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-5-target_S20_RL-merged.fastq.gz",
                   experiment = "G4", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-5-pegRNA_S5_RL-merged.fastq.gz",
                   experiment = "G4", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)
#45-4-6
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-6-target_S21_RL-merged.fastq.gz",
                   experiment = "G4", "ACTGAGCACG", "TGATGGCAGA", library_HEK3_target)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20211023_hiseq_walk-up_XXX_goose-mmr/45-4-6-pegRNA_S6_RL-merged.fastq.gz",
                   experiment = "G4", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGTT", library_HEK3)

#RT
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20220126_walkup_426/MiSeq_walk_up_436-322074858/FASTQ_Generation_2022-01-26_04_04_52Z-518947429/45_2_4_RT_L001-ds.0ca3026b04d6433c8eaf793de5f50b06/45-2-4-RT_S7_L001_R1_001.fastq.gz",
                   experiment = "RT", "TTTCCTCTGCCATCA", "CGTGCTCAGTCT", library_HEK3_RT)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20220126_walkup_426/45-2-5-RT_S8_L001_R1_001.fastq.gz",
                   experiment = "RT", "TTTCCTCTGCCATCA", "CGTGCTCAGTCT", library_HEK3_RT)
read_match_pegRNAs(workdir = "/Users/jk24/OneDrive/PhD/Prime_editing_efficiencies", datapath = "/raw_data/NGS/20220126_walkup_426/45-2-6-RT_S9_L001_R1_001.fastq.gz",
                   experiment = "RT", "TTTCCTCTGCCATCA", "CGTGCTCAGTCT", library_HEK3_RT)

# Canadian Goose target workdir = "/Volumes/drive/walkup_441"
# Canadian Goose target workdir = "/Volumes/drive/walkup_448"

# 18 nucleotide insertions
# HEK3_1
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-1-target_S33_L001_R1_001.fastq.gz",
                   experiment = "G5", "CAGACTGAGCACG", "TGATGGCA.AGGAAA", library_canadian_SP1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-2-target_S34_L001_R1_001.fastq.gz",
                   experiment = "G5", "CAGACTGAGCACG", "TGATGGCA.AGGAAA", library_canadian_SP1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-3-target_S35_L001_R1_001.fastq.gz",
                   experiment = "G5", "CAGACTGAGCACG", "TGATGGCA.AGGAAA", library_canadian_SP1_target)

# HEK3_2
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-4-target_S36_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-5-target_S37_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-6-target_S38_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2_target)
# HEK3_3
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-7-target_S39_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-8-target_S40_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-9-target_S41_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3_target)
# HEK3_4
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-10-target_S42_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-11-target_S43_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-12-target_S44_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4_target)
# HEK3_5
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-13-target_S67.fastq.gz",
                   experiment = "G5", "TGCCATGCCAGCT", "AAGTGGCTTTGGAGT", library_canadian_SP5_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-14-target_S68.fastq.gz",
                   experiment = "G5", "TGCCATGCCAGCT", "AAGTGGCTTTGGAGT", library_canadian_SP5_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-15-target_S69.fastq.gz",
                   experiment = "G5", "TGCCATGCCAGCT", "AAGTGGCTTTGGAGT", library_canadian_SP5_target)
# HEK3_6
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-16-target_S4_L001_R2_001_rc.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-17-target_S5_L001_R2_001_rc.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-18-target_S6_L001_R2_001_rc.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6_target)
# HEK3_6_2
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-16-target_S70.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCG.GGATTTGA", library_canadian_SP6_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-17-target_S71.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCG.GGATTTGA", library_canadian_SP6_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-18-target_S72.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCG.GGATTTGA", library_canadian_SP6_target)


# pegRNAs
# HEK3_1
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-1-pegRNA_S15_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTG", library_canadian_SP1)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-2-pegRNA_S16_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTG", library_canadian_SP1)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-3-pegRNA_S17_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTG", library_canadian_SP1)
# HEK3_2
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-4-pegRNA_S18_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-5-pegRNA_S19_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-6-pegRNA_S20_L001_R1_001.fastq.gz",
                   experiment = "G5", "AAACTGAGGCCAGAA", "AGTGATGGAGCTT", library_canadian_SP2)
# HEK3_3
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-7-pegRNA_S21_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-8-pegRNA_S22_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-9-pegRNA_S23_L001_R1_001.fastq.gz",
                   experiment = "G5", "TTGTTATGTCCTTTC", "ATCCTAGCAACTT", library_canadian_SP3)
# HEK3_4
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-10-pegRNA_S24_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-11-pegRNA_S25_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-12-pegRNA_S26_L001_R1_001.fastq.gz",
                   experiment = "G5", "TGGAAAACCCCAAAG", "GGTCCCTCTGACT", library_canadian_SP4)
# HEK3_5
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-13-pegRNA_S27_L001_R1_001.fastq.gz",
                   experiment = "G5", "ACTCCAAAGCCACTT", "AGCTGGCATGGCA", library_canadian_SP5)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-14-pegRNA_S28_L001_R1_001.fastq.gz",
                   experiment = "G5", "ACTCCAAAGCCACTT", "AGCTGGCATGGCA", library_canadian_SP5)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-15-pegRNA_S29_L001_R1_001.fastq.gz",
                   experiment = "G5", "ACTCCAAAGCCACTT", "AGCTGGCATGGCA", library_canadian_SP5)
# HEK3_6
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-16-pegRNA_S30_L001_R1_001.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-17-pegRNA_S31_L001_R1_001.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6)
read_match_pegRNAs(workdir = workdir, datapath = "/45-5-18-pegRNA_S32_L001_R1_001.fastq.gz",
                   experiment = "G5", "CATCTCCTGCCCAAA", "TGCGAGGATTTGA", library_canadian_SP6)

workdir = "/Volumes/drive/walkup_441"
workdir = "/Volumes/drive/walkup_XXX"
workdir = "/Volumes/drive/hiseq_walkup_180"
workdir = "/Users/jk24/Desktop/hiseq_1"

# HEK3_epegRNA
read_match_pegRNAs(workdir = workdir, datapath = "/45-10-1_pegRNA_S7_L001.assembled.fastq",
                   experiment = "G6", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGCG", library_epeg)
read_match_pegRNAs(workdir = workdir, datapath = "/45-10-2_pegRNA_S8_L001.assembled.fastq",
                   experiment = "G6", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGCG", library_epeg)
read_match_pegRNAs(workdir = workdir, datapath = "/45-10-3_pegRNA_S9_L001.assembled.fastq",
                   experiment = "G6", "TTTCCTCTGCCATCA", "CGTGCTCAGTCTGCG", library_epeg)

read_match_pegRNAs(workdir = workdir, datapath = "/45-10-1_target_S16_L001.assembled.fastq",
                   experiment = "G6", "ACTGAGCACG", "TGATGGCAGA", library_epeg_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-10-2_target_S17_L001.assembled.fastq",
                   experiment = "G6", "ACTGAGCACG", "TGATGGCAGA", library_epeg_target)
read_match_pegRNAs(workdir = workdir, datapath = "/45-10-3_target_S18_L001.assembled.fastq",
                   experiment = "G6", "ACTGAGCACG", "TGATGGCAGA", library_epeg_target)


# ==== BRENT (Set 2 sequences) ====
# Brent Library
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-CLYBL-a_S1_L001_R1_001.fastq.gz",
                   experiment = "Library", "AGGTCCTTAG", "TGACAGTCTG", filter(library_brent, subpool == "small_CLYBL"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-CLYBL-b_S2_L001_R1_001.fastq.gz",
                   experiment = "Library", "AGGTCCTTAG", "TGACAGTCTG", filter(library_brent, subpool == "small_CLYBL"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-EMX1-a_S3_L001_R1_001.fastq.gz",
                   experiment = "Library", "GGAGCCCTTC", "TTCTTCTGCT", filter(library_brent, subpool == "small_EMX1"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-EMX1-b_S4_L001_R1_001.fastq.gz",
                   experiment = "Library", "GGAGCCCTTC", "TTCTTCTGCT", filter(library_brent, subpool == "small_EMX1"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-FANCF-a_S5_L001_R1_001.fastq.gz",
                   experiment = "Library", "CGATCCAGGT", "GCTGCAGAAG", filter(library_brent, subpool == "small_FANCF"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-FANCF-b_S6_L001_R1_001.fastq.gz",
                   experiment = "Library", "CGATCCAGGT", "GCTGCAGAAG", filter(library_brent, subpool == "small_FANCF"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-HEK3-a_S7_L001_R1_001.fastq.gz",
                   experiment = "Library", "TCTGCCATCA", "CGTGCTCAGT", filter(library_brent, subpool == "small_HEK3"))
read_match_pegRNAs(workdir = workdir, datapath = "/Brent-HEK3-b_S8_L001_R1_001.fastq.gz",
                   experiment = "Library", "TCTGCCATCA", "CGTGCTCAGT", filter(library_brent, subpool == "small_HEK3"))


# Brent pegRNAs
# L1
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-7-pegRNA_S4_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-8-pegRNA_S5_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-9-pegRNA_S6_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-10-pegRNA_S7_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-11-pegRNA_S8_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-12-pegRNA_S9_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-13-pegRNA_S10_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-14-pegRNA_S11_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-15-pegRNA_S12_L001.assembled.fastq.gz", experiment = "Brent", library_brent)
# L2
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-7-pegRNA_S4_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-8-pegRNA_S5_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-9-pegRNA_S6_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-10-pegRNA_S7_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-11-pegRNA_S8_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-12-pegRNA_S9_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-13-pegRNA_S10_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-14-pegRNA_S11_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)
read_match_pegRNAs2(workdir = workdir, datapath = "/pegRNA/45-6-15-pegRNA_S12_L002.assembled.fastq.gz", experiment = "Brent_L2", library_brent)

# Brent targets
# FANCF
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-7-FANCF_S33.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-8-FANCF_S37.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-9-FANCF_S41.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-10-FANCF_S45.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-11-FANCF_S49.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-12-FANCF_S53.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-13-FANCF_S57.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-14-FANCF_S61.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
read_match_pegRNAs(workdir = workdir, datapath = "/FANCF/45-6-15-FANCF_S65.fastq.gz",
                   experiment = "Brent", "CTTCTGCAGC", "ACCTGGATCG", library_brent_FANCF_target)
# CLYBL
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-7-CLYBL_S31.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-8-CLYBL_S35.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-9-CLYBL_S39.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-10-CLYBL_S43.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-11-CLYBL_S47.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-12-CLYBL_S51.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-13-CLYBL_S55.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-14-CLYBL_S59.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
read_match_pegRNAs(workdir = workdir, datapath = "/CLYBL/45-6-15-CLYBL_S63.fastq.gz",
                   experiment = "Brent", "AGGTCCTTAG", "TGACAGTCTG", library_brent_CLYBL_target)
# EMX1
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-7-EMX1_S32.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-8-EMX1_S36.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-9-EMX1_S40.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-10-EMX1_S44.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-11-EMX1_S48.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-12-EMX1_S52.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-13-EMX1_S56.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-14-EMX1_S60.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)
read_match_pegRNAs(workdir = workdir, datapath = "/EMX1/45-6-15-EMX1_S64.fastq.gz",
                   experiment = "Brent", "AGCAGAAGAA", "GAAGGGCTCC", library_brent_EMX1_target)

# HEK3 L1
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-7-HEK3_S34_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-8-HEK3_S38_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-9-HEK3_S42_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-10-HEK3_S46_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-11-HEK3_S50_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-12-HEK3_S54_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-13-HEK3_S58_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-14-HEK3_S62_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-15-HEK3_S66_L002.assembled.fastq.gz",
                   experiment = "Brent", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
# HEK3 L2
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-7-HEK3_S34_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-8-HEK3_S38_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-9-HEK3_S42_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCA.A", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-10-HEK3_S46_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-11-HEK3_S50_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-12-HEK3_S54_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-13-HEK3_S58_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-14-HEK3_S62_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)
read_match_pegRNAs(workdir = workdir, datapath = "/HEK3/45-6-15-HEK3_S66_L002.assembled.fastq.gz",
                   experiment = "Brent_L2", "ACTGAGCACG", "TGATGGCAGA", library_brent_HEK3_target)

# barnacle workdir = "/Volumes/drive/walkup_180_hiseq"
# ACTB
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-ACTB_S13.fastq.gz",
                   experiment = "Barnacle", "GCTCACCATG", "GATGATGATA", filter(library_barnacle_target, target == "ACTB" | target == "ACTB_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-ACTB_S19.fastq.gz",
                   experiment = "Barnacle", "GCTCACCATG", "GATGATGATA", filter(library_barnacle_target, target == "ACTB" | target == "ACTB_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-ACTB_S25.fastq.gz",
                   experiment = "Barnacle", "GCTCACCATG", "GATGATGATA", filter(library_barnacle_target, target == "ACTB" | target == "ACTB_rc"))
# HEK3
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-HEK3_S14.fastq.gz",
                   experiment = "Barnacle", "ACTGAGCACG", "TGATGGCA.A", filter(library_barnacle_target, target == "HEK3"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-HEK3_S20.fastq.gz",
                   experiment = "Barnacle", "ACTGAGCACG", "TGATGGCA.A", filter(library_barnacle_target, target == "HEK3"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-HEK3_S26.fastq.gz",
                   experiment = "Barnacle", "ACTGAGCACG", "TGATGGCA.A", filter(library_barnacle_target, target == "HEK3"))
# LMNB1
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-LMNB1_S15.fastq.gz",
                   experiment = "Barnacle", "GCCCGCCATG", "GCGACTGCGA", filter(library_barnacle_target, target == "LMNB1" | target == "LMNB1_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-LMNB1_S21.fastq.gz",
                   experiment = "Barnacle", "GCCCGCCATG", "GCGACTGCGA", filter(library_barnacle_target, target == "LMNB1" | target == "LMNB1_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-LMNB1_S27.fastq.gz",
                   experiment = "Barnacle", "GCCCGCCATG", "GCGACTGCGA", filter(library_barnacle_target, target == "LMNB1" | target == "LMNB1_rc"))
# NOLC1
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-NOLC1_S16.fastq.gz",
                   experiment = "Barnacle", "CTGGAGGATG", "GCGGACGCCG", filter(library_barnacle_target, target == "NOLC1"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-NOLC1_S22.fastq.gz",
                   experiment = "Barnacle", "CTGGAGGATG", "GCGGACGCCG", filter(library_barnacle_target, target == "NOLC1"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-NOLC1_S28.fastq.gz",
                   experiment = "Barnacle", "CTGGAGGATG", "GCGGACGCCG", filter(library_barnacle_target, target == "NOLC1"))
# RNF1
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-RNF2_S17.fastq.gz",
                   experiment = "Barnacle", "TAGTCATTAC", "CTGAGGTGTT", filter(library_barnacle_target, target == "RNF1"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-RNF2_S23.fastq.gz",
                   experiment = "Barnacle", "TAGTCATTAC", "CTGAGGTGTT", filter(library_barnacle_target, target == "RNF1"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-RNF2_S29.fastq.gz",
                   experiment = "Barnacle", "TAGTCATTAC", "CTGAGGTGTT", filter(library_barnacle_target, target == "RNF1"))
# TP53
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-1-TP53_S18.fastq.gz",
                   experiment = "Barnacle", "CACTGCCATG", "GAGGAGCCGC", filter(library_barnacle_target, target == "TP53" | target == "TP53_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-2-TP53_S24.fastq.gz",
                   experiment = "Barnacle", "CACTGCCATG", "GAGGAGCCGC", filter(library_barnacle_target, target == "TP53" | target == "TP53_rc"))
read_match_pegRNAs(workdir = workdir, datapath = "/barnacle/45-6-3-TP53_S30.fastq.gz",
                   experiment = "Barnacle", "CACTGCCATG", "GAGGAGCCGC", filter(library_barnacle_target, target == "TP53" | target == "TP53_rc"))


# barnacle pegRNA
read_match_pegRNAs2(workdir = workdir, datapath = "/barnacle/45-6-1-pegRNA_S1.fastq.gz",
                    experiment = "Barnacle", library_barnacle)
read_match_pegRNAs2(workdir = workdir, datapath = "/barnacle/45-6-2-pegRNA_S2.fastq.gz",
                    experiment = "Barnacle", library_barnacle)
read_match_pegRNAs2(workdir = workdir, datapath = "/barnacle/45-6-3-pegRNA_S3.fastq.gz",
                    experiment = "Barnacle", library_barnacle)

read_match_pegRNAs2(workdir = workdir, datapath = "/Barnacle-Goose-a_S9_L001_R1_001.fastq.gz",
                    experiment = "Library", library_barnacle)
read_match_pegRNAs2(workdir = workdir, datapath = "/Barnacle-Goose-b_S10_L001_R1_001.fastq.gz",
                    experiment = "Library", library_barnacle)





