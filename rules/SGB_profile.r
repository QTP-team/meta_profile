library(tidyverse)
library(readxl)

args=commandArgs(T)

infile <-  args[1]
outfile <- args[2]
genomeInfo_file <- args[3]
id_cut <- args[4]
reads_cut <- args[5]

genomeInfo <- read_tsv(genomeInfo_file)

SGB_size <- genomeInfo %>%   
  group_by(SGB_ID) %>%
  summarise(SGB_size = sum(Contig_length))

sel_alignment <- function(infile, id_cutoff = 0.95) {
  #infile = "0.data/mock.10M.map2SGB.aligned.reads.tsv"
  in.dset <- read_tsv(infile)
  uniq.reads <- in.dset %>%  
    group_by(Read_ID) %>%
    summarise(Read_n = n()) %>%
    filter(Read_n == 1) 
        
  uniq.filtered.dset <- in.dset %>%
    semi_join(uniq.reads, by = "Read_ID") %>%   
    filter(Identity >= id_cutoff)
  
  filter.out.n = nrow(uniq.reads)-nrow(uniq.filtered.dset)
  #print(paste0("filter out ",filter.out.n, " reads which identity lower than ", id_cutoff))
  return(uniq.filtered.dset)
}

profile <- function(uniq.filtered.dset, reads_cutoff = 0){
  contig_r <- uniq.filtered.dset %>%
    group_by(Contig_ID) %>%
    summarise(Contig_reads = n()) %>%
    ungroup()
  
  SGB_r <- contig_r %>%
  left_join(select(genomeInfo, SGB_ID, Contig_ID), by = "Contig_ID") %>%
  group_by(SGB_ID) %>%
  summarise(SGB_reads = sum(Contig_reads)) %>%
  ungroup()
  
  SGB_r_filter <- SGB_r %>%
  filter(SGB_reads > reads_cutoff) %>%
  left_join(SGB_size, by = "SGB_ID") %>%
  mutate(`Fragment/gSize` = SGB_reads/SGB_size,
         FPKM = (1e3 * 1e6 * SGB_reads) / (sum(SGB_reads) * SGB_size),
         rel_ab = `Fragment/gSize` / sum(`Fragment/gSize`) * 100)
  
  #print(paste0(nrow(SGB_r), " SGBs are found in the sample, ", nrow(SGB_r_filter), " SGBs which alignments higher than ", reads_cutoff, " are retained."))
  return(SGB_r_filter)
}

uniq.filtered.dset <- sel_alignment(infile, as.numeric(id_cut))
SGB_profile <- profile(uniq.filtered.dset, as.numeric(reads_cut))

write_tsv(SGB_profile, outfile)