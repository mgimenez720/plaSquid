#!/usr/bin/env Rscript

  args = commandArgs(trailingOnly=TRUE)

   library(readr)
   library(tidyverse)
   library(Biostrings)

  mbt = args[1] #Minidist_result.tsv
  rbt = args[2] #Plasmid_Report.tsv
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
  ltb <-  inner_join (tbm, tns, by = c("Contig" = "names"))


  ftb <- full_join(ltb, rtb, by = c("Contig"))

  ftb1 <- ftb[,c("names","Contig","Sim_dist","plsdb_match","Match_length","Rep_type", "MOB_group", "RIP_domain", "Contig_length.y")]
  
  ftb2 <- full_join(ftb1, tns, by= c("Contig" = "names"))
  
  ftb3 <- ftb2[,c("Contig.y","Contig","Sim_dist","plsdb_match","Match_length","Rep_type", "MOB_group", "RIP_domain", "Contig_length.y")]
  
  names(ftb3) <- c("Contig", "name", "Sim_dist","plsdb_match","Match_length","Inc_group", "MOB_group", "RIP_domain", "Contig_length")


  ftb4 <- subset.data.frame(ftb3, subset = ftb3$Sim_dist>45 | ftb3$Inc_group != "NA" | ftb3$MOB_group != "NA" | ftb3$RIP_domain != "NA")
  
   hts <- ftb4$name
   idx <- which(nfs %in% hts)
   fst1 <- fst[idx]

   writeXStringSet(fst1, "Result.fasta")
   write_delim(ftb4,"Result.tsv", delim= "\t")


