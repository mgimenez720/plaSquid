#!/bin/bash/env Rscript

 args = commandArgs(trailingOnly=TRUE)

 tbr = args[1]
 tbm = args[2]
 asm = args[3]
 
 library(tidyverse)
 library(Biostrings)
 
 #Load data
 
 rps <- read_delim(tbr, delim="\t")
 mbs <- read_delim(tbm, delim="\t")
 mtg <- readDNAStringSet(asm)
 
 #Collect contigs
 
 rpc <- rps$contig
 mbc <- mbs$contig
 
 cnts <- unique(c(rpc, mbc))
 
 #Extract plasmidic contigs

 nms <- names(mtg)
 idx <- which(nms %in% cnts)
 
 hit <- mtg[idx]
 
 
 #Parsing output tables
 
 names(rps) <- c("Rep_ORF", "Rep_type", "Rep_score", "Rep_len", "contig")
 names(mbs) <- c("Mob_ORF", "Mob_len", "MOB_group", "MOB_score", "alifrom", "alito", "contig")
 
 rps1 <- rps %>% 
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 mbs1 <- mbs %>%
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 ftb <- full_join(rps1, mbs1, by = "contig")
 
 fct <- ftb$contig
 nht <- names(hit)
 len <- numeric(0)
 
 for (i in 1:length(fct)){
   
   ct <- fct[i]
   idx <- which(nht == ct)
   cnl <- width(hit[idx])
   
   len <- c(len, cnl)
   
 }
 
 ftb$contig_length <- len
 
 ftb1 <- ftb[,c(1,3,8,12,2,6)]
 
 #Writing final results
 
 writeXStringSet(hit, "Plasmids_contigs.fasta")
 write_delim(ftb1, "Plasmid_Report.tsv")
 
 