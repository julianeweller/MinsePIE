# Loading libraries
library(ShortRead)
library(spgs)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
# Arguments to be supplied
# 1 library path (a file with all library insert sequences)
# 2 data path (a path to the sequencing fastq file)
# 3 outdir (directory to write the results into)
library_mut <- read_csv(args[1])
library_pattern <- paste(filter(library_mut, nchar(sequence_original) > 2)$sequence_original, collapse = "|")

read_sequences <- function(datapath, outdir) {
  
  # Extract the useful parts of the file name
  filename = sub("_S.*$", "", basename(datapath))
  print(filename)
  
  print("Reading fastq files")
  reads <- readFastq(datapath)
  reads <- as.character(reads@sread)
  # Trim reads
  reads_trim <- str_extract(reads, "GTCCGAGCAG(.*?)TAGGGTGGGC")
  
  # Trim reads
  print("Counting and matching pegRNAs")
  count_table <-   count_table <- tibble(sequence = reads_trim) %>%
    filter(nchar(sequence) < 60) %>% 
    # First collapse reads
    group_by(sequence) %>%
    summarize(count = n()) %>%
    # Identifying pegRNAs that were not detected in sequencing data
    mutate(count = ifelse(is.na(count), 0, count), sequence_trim = substr(sequence, 16, nchar(sequence)-15), sample = filename) %>%
    select(sequence, count, sequence_trim, sample) %>%
    filter(!is.na(sequence))

  print("Detecting sequences with library inserts")
  # Some clean up. Fixing duplicatred sequences in the library. Normalization
  count_table <- count_table %>% 
    mutate(sample = filename, length = nchar(sequence),
           library = str_detect(sequence_trim, library_pattern),
           trim_in_lib = sequence_trim %in% library_mut$sequence_original) %>%
    mutate(library = ifelse(trim_in_lib == F, F, library)) %>%
    arrange(desc(count))
  
  write_tsv(count_table, paste0(outdir, filename, "_nicking_sequences.tsv"))
  write_tsv(tibble("total_reads" = length(reads), "deletions" = sum(count_table$count), "sample" = filename), paste0(outdir, filename, "_nicking_metrics.tsv"))
  

}

read_sequences(datapath = args[2], outdir = args[3])