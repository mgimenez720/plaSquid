#!/usr/bin/env Rscript

library(Biostrings)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

tbb = args[1]
dnn = args[2]

dna <- readDNAStringSet(dnn)

tab <- read_delim(tbb, delim = "\t")

cnt <- tab$Contig
idx <- numeric()

for (i in 1:length(cnt)) {
  
  cn <- cnt[i]
  id <- strsplit(cn, "-")[[1]][2]
  
  paste0()
    
  idx <- c(idx,id)
  
  
}
  
idx <- as.numeric(idx)
dnr <- dna[idx]
names(dnr)<-tab$Contig
writeXStringSet(dnr, "Result.fasta")

colnames(tab) <- c("Contig", "RIP_domain", "MOB_group", "Rep_type", "contig_length")
write_delim(tab, "Result.tsv", delim = "\t")
