#!/bin/bash/env Rscript

args = commandArgs(trailingOnly=TRUE)

tab = args[1]         # Input paf table
fst = args[2]         # Input fasta file

 library(readr)
 library(Biostrings)
 library(dplyr)
 
 tab1 <- read_delim(tab, delim = "\t", col_names = FALSE)
 colnames(tab1) <- c("Contigs", "Sim-dist", "plsdb_match", "match_length")
 tab2 <- arrange(tab1, Contigs)
 
 idx <- which(tab2$"Sim-dist" >= 45)
 
 tab3 <- tab2[idx,]
 plsn <- tab3$Contigs
 
 #Read fasta and change names to Contigs
 
 fas <- readDNAStringSet(filepath = fst)
 nff <- names(fas)
 
 
 nms <- character(0)
 
 for (i in 1:length(fas)) {
   
   nm <- paste0("Contig-",i)
   
   nms <- c(nms, nm) 
   
 }
 
   names(fas) <- nms 
   id1 <- which(nms %in% plsn)

   fsf <- fas[id1]
   nfp <- nff[id1]    
   lfp <- width(fsf)
   tab3$Contigs <- nfp
   tab3$Contig_length <- lfp
   
   names(fsf) <- nfp
   
   writeXStringSet(fsf, filepath = "Minidist_contigs.fasta")
   write_delim(tab3, path = "Minidist_result.tsv", delim = "\t")
   
   