#!/usr/bin/env Rscript

  args = commandArgs(trailingOnly=TRUE)

   library(readr)
   library(tidyverse)
   library(Biostrings)

  mbt = args[1] #Minidist_result.tsv
  fsn = args[2] #test.fasta

   
   tbm <- read_delim(mbt, delim ="\t")

   fst <- readDNAStringSet(fsn)
   nfs <- names(fst)


    names(tbm) <- c("name","Sim_dist","plsdb_match","Match_length","Contig_length")
    
    tbm1 <- subset.data.frame(tbm, subset = tbm$Sim_dist > 45)
  
    hts <- tbm1$name
    idx <- which(nfs %in% hts)
    fst1 <- fst[idx]

    writeXStringSet(fst1, "Result.fasta")
    write_delim(tbm1,"Result.tsv", delim= "\t")


#Falta agregar un proceso con este script cuando llamo solo a Minidist y agregar una carpeta de output! como hace sum info???