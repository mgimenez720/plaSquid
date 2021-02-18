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

  rtb <-  full_join(tns, tbr, by = "Contig")
  names(rtb) <- c("Contig", "names", "Rep_type", "MOB_group", "RIP_domain", "Contig_length")
  names(tbm) <- c("Contig","Sim_dist","plsdb_match","Match_length","Contig_length")

  ftb <- full_join(tbm, rtb, by = c("Contig", "Contig_length"))

  ftb1 <- ftb[,c(1,2,3,4,5,7,8,9)]


   hts <- ftb1$Contig
   idx <- which(hts %in% nfs)
   fst1 <- fst[idx]

   writeXStringSet(fst1, "Result.fasta")
   write_delim(ftb1,"Result.tsv", delim= "\t")


