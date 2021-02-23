#!/usr/bin/env Rscript

  args = commandArgs(trailingOnly=TRUE)

   library(readr)
   library(tidyverse)
   library(Biostrings)

  mbt = args[1] #Minidist_result.tsv
  rbt = args[2] #Plasmid_report.tsv
  fsn = args[3] #test.fasta

   
   tbm <- read_delim(mbt, delim ="\t")
   tbr <- read_delim(rbt, delim = "\t")

   fst <- readDNAStringSet(fsn)
   nfs <- names(fst)

   cnts <- character(0)
   for (i in 1:length(nfs)){

     cs <- paste0("contig-",i)

     cnts <- c(cnts, cs)

   }

     tns <- tibble(names = nfs,
                   Contig = cnts)

  rtb <-  inner_join(tns, tbr, by = "Contig")
  names(rtb) <- c("Contig", "names", "Rep_type", "MOB_group", "RIP_domain", "Contig_length")
  
  
  names(tbm) <- c("Contig","Sim_dist","plsdb_match","Match_length","Contig_length")
  ltb <-  inner_join(tns, tbm, by = c("names" = "Contig"))


  ftb <- full_join(ltb, rtb, by = c("Contig"))

  ftb1 <- ftb[,c("names.x","Contig","Sim_dist","plsdb_match","Match_length","Inc_group", "MOB_group", "RIP_domain", "contig_length")]
  names(ftb1) <- c("Contig", "name", "Sim_dist","plsdb_match","Match_length","Inc_group", "MOB_group", "RIP_domain", "Contig_length")


   hts <- ftb1$Contig
   idx <- which(hts %in% nfs)
   fst1 <- fst[idx]

   writeXStringSet(fst1, "Result.fasta")
   write_delim(ftb1,"Result.tsv", delim= "\t")


