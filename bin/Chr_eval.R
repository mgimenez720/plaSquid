#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


res = args[2]
dna = args[1]

library(tidyverse)
library(Biostrings)

    hmt <- as_tibble(read.table(res, header=TRUE, sep = "\t", comment.char = '#', 
                                 col.names= c("Contig", "name", "Sim_dist", "plsdb_match","Match_length",
                                             "RIP_domain", "MOB_group", "Rep_type", "Contig_length"),
                                fill = TRUE, stringsAsFactors = FALSE))

    

    idx <- which(is.na(hmt$plsdb_match))

    dnn <- readDNAStringSet(dna)
    names(dnn) <- hmt$Contig
    
    dn1 <- dnn[idx]

    writeXStringSet(dn1, "Chr_eval.fasta")
   
      
   