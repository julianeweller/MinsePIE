# 2021-11-08
# This script takes fastq sequencing files as input and generates tables with SNP frequencies and a list of all sequences, their frequencies, and whether or not they correspond to a library insertion as output

# Loading libraries
library(ShortRead)
library(spgs)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
# Arguments to be supplied
# 1 library path (a file with all library insert sequences)
# 2 data path (a path to the sequencing fastq file)
# 3 experiment (name of the experiment)
# 4 leftcons (a string of the sequence the script should match at the 5' end of the amplicon)
# 5 rightcons (a string of the sequence the script should match at the 3' end of the amplicon)
# 6 outdir (directory to write the results into)
# 7 preserved (a string of the wild type target site sequence)

# Importing files
library_mut <- read_csv(args[1])

# Defining a set of functions
rc <- function(x) {toupper(spgs::reverseComplement(x))}

count_nucleotides <- function(x, outdir, filename) {
  substrings = matrix(nrow = nrow(x), ncol = 50)
  substringsr = matrix(nrow = nrow(x), ncol = 50)
  
  for(i in 1:50) {
    substrings[, i] = substring(x$sequence, i+10, i+10)
  }
  for(i in 1:50) {
    substringsr[, i] = substring(x$sequence, 
                                 nchar(x$sequence)-i-9, nchar(x$sequence)-i-9)
  }
  
  substringsA = colSums((substrings == "A")*x$count)
  substringsC = colSums((substrings == "C")*x$count)
  substringsG = colSums((substrings == "G")*x$count)
  substringsT = colSums((substrings == "T")*x$count)
  substringsAr = colSums((substringsr == "A")*x$count)
  substringsCr = colSums((substringsr == "C")*x$count)
  substringsGr = colSums((substringsr == "G")*x$count)
  substringsTr = colSums((substringsr == "T")*x$count)
  
  nuc_frequency = tibble("position" = c(-49:0), "A" = substringsA, "C" = substringsC, "G" = substringsG, "T" = substringsT) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "count")
  nuc_frequencyr = tibble("position" = seq(50, 1, by = -1), "A" = substringsAr, "C" = substringsCr, "G" = substringsGr, "T" = substringsTr) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "count")
  nuc_frequencyc = bind_rows(nuc_frequency, nuc_frequencyr) %>% group_by(position) %>%
    mutate(count = count/sum(count, na.rm = T)*100, sample = filename)
  write_tsv(nuc_frequencyc, paste0(outdir, filename, "_counts.tsv"))
}

find_indels <- function(x, outdir, filename, preserved) {
  indel_table <- filter(x, library == F) %>%
    mutate(preserved_target = str_detect(sequence, preserved)) %>%
    filter(preserved_target == F & !str_detect(sequence, "N")) %>%
    filter(length != 120) %>%
    filter(!nchar(sequence_trim) %in% c(1,2))
  write_tsv(indel_table, paste0(outdir, filename, "_indel_sequences.tsv"))
  indel_count <- sum(indel_table$count)/sum(x$count)*100
  write_tsv(tibble("percIndel" = indel_count, "sample" = filename), paste0(outdir, filename, "_indel_counts.tsv"))
}

read_sequences <- function(datapath, experiment, leftcons, rightcons, outdir) {
  
  # Extract the useful parts of the file name
  filename = sub("_S.*$", "", basename(datapath))
  print(filename)
  
  print("Reading fastq files")
  reads <- readFastq(datapath)
  reads <- as.character(reads@sread)
  # Trim reads
  reads_trim <- str_extract(reads, paste0(leftcons, "(.*?)", rightcons))
  
  print("Counting and matching pegRNAs")
  count_table <- tibble(sequence = reads_trim) %>%
    # First collapse reads
    group_by(sequence) %>%
    summarize(count = n()) %>%
    # Identifying pegRNAs that were not detected in sequencing data
    mutate(count = ifelse(is.na(count), 0, count), sequence_trim = substr(sequence, 61, nchar(sequence)-60), sample = filename) %>%
    select(sequence, count, sequence_trim, sample) %>%
    filter(!is.na(sequence))
  
  # Remove 1 and 2 nt sequences from library to prevent faulty matching
  library_pattern <- paste(filter(library_mut, nchar(sequence_original) > 2)$sequence_original, collapse = "|")
  
  print("Detecting sequences with library inserts")
  # Some clean up. Fixing duplicatred sequences in the library. Normalization
  count_table <- count_table %>% 
    mutate(sample = filename, length = nchar(sequence),
           library = str_detect(sequence_trim, library_pattern),
           trim_in_lib = sequence_trim %in% library_mut$sequence_original) %>%
    mutate(library = ifelse(trim_in_lib == F, F, library)) %>%
    arrange(desc(count))
  
  write_tsv(count_table, paste0(outdir, filename, "_sequences.tsv"))
  
  print("Detecting indels")
  find_indels(count_table, outdir = args[6], filename = filename, preserved = args[7])
  
  print("Count nucleotides")
  count_nucleotides(count_table, filename = filename, outdir = args[6])
  count_nucleotides(filter(count_table, library == T), filename = paste0(filename, "_library"), outdir = args[6])

}

# Running functions
# Example to run in a unix command line environment: Rscript --vanilla ./scripts/find_SNVs.R "./files/library.csv" "./fastq/45-2-1-target_S19_LRmerged.fastq.gz" "M" "GTCTCCAAGG" "AGCACCTGGG" "./prc_data/mutation_data/" "TCCCTTCTGCAGCACCTGGATCGCTTTTCCGAGCTTCTGGCGGTCTC"

read_sequences(datapath = args[2], experiment = args[3], leftcons = args[4], rightcons = args[5], outdir = args[6])