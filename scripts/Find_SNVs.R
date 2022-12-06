# 2022-11-11
# This script takes fastq sequencing files as input and generates tables with SNV frequencies and a list of all undesired sequences, their frequencies, type, and whether or not they correspond to a library insertion as output

# Loading libraries
library(ShortRead)
library(spgs)
library(tidyverse)
library(fuzzyjoin)

args = commandArgs(trailingOnly=TRUE)
# Arguments to be supplied
# 1 library path (a file with all library insert sequences)
# 2 data path (a path to the sequencing fastq file)
# 3 experiment (name of the experiment)
# 4 leftcons (a string of the sequence the script should match at the 5' end of the amplicon)
# 5 rightcons (a string of the sequence the script should match at the 3' end of the amplicon)
# 6 outdir (directory to write the results into)
# 7 preserved (a string of the wild type target site sequence, extending 3 nt from HA and RTT)
# 8 HA (a string of the HA)
# Example to run in a unix command line environment: 
# Rscript --vanilla find_SNVs.R "./files/library.csv" "./fastq/sample.fastq" "mutation" "CCTTGGGGCC" "AGCTTTTCCT" "./mut_data/" "GCCCAGACTGAGCACGTGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCAGAG" "TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA"


# Importing files
library_mut <- read_csv(args[1])

# Defining a set of functions
rc <- function(x) {toupper(spgs::reverseComplement(x))}

count_nucleotides <- function(x, outdir, filename) {
  # This function computes the fraction of nucleotides at each position along the sequence, going inwards from the outside and ending at the nick site
  substrings = matrix(nrow = nrow(x), ncol = 13)
  substringsr = matrix(nrow = nrow(x), ncol = 50)
  
  for(i in 1:13) {
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
  
  nuc_frequency = tibble("position" = c(-12:0), "A" = substringsA, "C" = substringsC, "G" = substringsG, "T" = substringsT) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "count")
  nuc_frequencyr = tibble("position" = seq(50, 1, by = -1), "A" = substringsAr, "C" = substringsCr, "G" = substringsGr, "T" = substringsTr) %>%
    pivot_longer(cols = c("A", "C", "G", "T"), names_to = "nucleotide", values_to = "count")
  nuc_frequencyc = bind_rows(nuc_frequency, nuc_frequencyr) %>% group_by(position) %>%
    mutate(count = count/sum(count, na.rm = T)*100, sample = filename)
  write_tsv(nuc_frequencyc, paste0(outdir, filename, "_counts.tsv"))
}

find_indels <- function(x, outdir, filename, preserved) {
  # This function applies multiple filters to the set of sequences to identify whether they are an undesired editing product
  indel_table <- filter(x, library == F) %>% # 1: filter out sequences with perfect library insertions
    mutate(preserved_target = str_detect(sequence, preserved)) %>% # 2: filter out sequences, that have a perfectly preserved wild type target site sequence, extending 3 nt from HA and RTT  
    filter(preserved_target == F & !str_detect(sequence, "N")) %>% # 3: filter out sequences with ambiguous sequencing bases (N)
    filter(length != 83) # 4: filter out sequences that have the length corresponding to a wild type sequence (these will be analyzed in the count_nucleotides function instead)
}

annotate_indels <- function(x, outdir, filename, HA, preserved) {
  # This function annotateds indels
  x_filtered <- filter(x, nchar(sequence_trim) > 10) %>% group_by(sequence_trim) %>% summarise() # approximately match library sequences > 10 nt to detect mutated insertions
  almost_matching <- stringdist_inner_join(library_mut, x_filtered, by = c("sequence_original" = "sequence_trim"), method = 'osa', max_dist = 3) # up to 3 mismatches are allowed
  
  # create a vector to catch deletions at the target site
  target_deletions <- vector(mode = "character")
  HA_deletions <- vector(mode = "character")
  PBS_deletions <- vector(mode = "character")
  
  for(i in 1:10) {
    target_deletions[i] = paste0(substr(preserved, 4, 16-i), substr(HA, i+1, nchar(HA)))
    HA_deletions[i] = paste0(substr(preserved, 4, 16), substr(HA, i+1, nchar(HA)))
    PBS_deletions[i] = paste0(substr(preserved, 4, 16-i), HA)
  }
  
  deletions <- paste(c(target_deletions, HA_deletions, PBS_deletions), collapse = "|")
  
  annotated_indels <- x %>% 
    left_join(almost_matching, by = "sequence_trim") %>%
    mutate(
      # First check if we see scaffold insertion
      outcome_type = 
        ifelse(str_detect(sequence, paste0(HA, "GCACC")), "scaffold_integration",
             # Then check if we see an approximate match (2 nt) to a library sequence
             ifelse(!is.na(sequence_original), "approximate_match",
                   # Check for presence of a tandem duplication involving the HA
                   ifelse(str_detect(sequence_trim, HA), "duplication",
                          # Check for deletion
                          ifelse(length < 83, "deletion", "other")))),
      outcome_type = ifelse(outcome_type == "deletion", ifelse(str_detect(sequence, deletions), "deletion_at_target", "deletion"), outcome_type)
    )
}

read_sequences <- function(datapath, experiment, leftcons, rightcons, outdir) {
  # A wrapper function that reads in the fastq file and runs individual functions on it
  
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
    mutate(count = ifelse(is.na(count), 0, count), sequence_trim = substr(sequence, 24, nchar(sequence)-60), sample = filename) %>%
    select(sequence, count, sequence_trim, sample) %>%
    filter(!is.na(sequence))
  
  # Remove 1 and 2 nt sequences from library to prevent faulty matching
  library_pattern <- paste(filter(library_mut, nchar(sequence_original) > 2)$sequence_original, collapse = "|")
  
  print("Detecting sequences with library inserts")
  # Some clean up. Fixing duplicated sequences in the library. Normalization
  count_table <- count_table %>% 
    mutate(sample = filename, length = nchar(sequence),
           library = str_detect(sequence_trim, library_pattern),
           trim_in_lib = sequence_trim %in% library_mut$sequence_original) %>%
    mutate(library = ifelse(trim_in_lib == F, F, library)) %>%
    arrange(desc(count))
  
  write_tsv(count_table, paste0(outdir, filename, "_sequences.tsv"))
  
  print("Detecting indels")
  indel_table <- find_indels(count_table, outdir = args[6], filename = filename, preserved = args[7])
  
  print("Annotate indels")
  annotated_indels <- annotate_indels(indel_table, outdir = args[6], filename = filename, preserved = args[7], HA = args[8])
  write_tsv(annotated_indels, paste0(outdir, filename, "_indel_sequences.tsv"))
  
  print("Summarising")
  summarised_indels <- annotated_indels %>% group_by(outcome_type) %>% 
    summarise(count = sum(count)) %>% 
    mutate(count_relative = count/sum(count_table$count)*100, "sample" = filename)
  write_tsv(summarised_indels, paste0(outdir, filename, "_indel_counts.tsv"))
  
  print("Count nucleotides")
  count_nucleotides(filter(count_table, length == 83 & library == F), filename = filename, outdir = args[6]) # only sequences that did not insert a sequence and show no indels
  count_nucleotides(filter(count_table, library == T), filename = paste0(filename, "_library"), outdir = args[6]) # only sequences in the library

}

# Running functions
read_sequences(datapath = args[2], experiment = args[3], leftcons = args[4], rightcons = args[5], outdir = args[6])
