#!/bin/bash/env Rscript

 args = commandArgs(trailingOnly=TRUE)

 tbr = args[1]
 tbd = args[2]
 tbm = args[3]
 asm = args[4]
 
 library(tidyverse)
 library(Biostrings)
 
 #Load data
 
 rps <- distinct(read_delim(tbr, delim="\t"), .keep_all = TRUE)
 rpd <- na.omit(read_delim(tbd, delim = "\t"))
 mbs <- read_delim(tbm, delim="\t")
 mtg <- readDNAStringSet(asm)
 
 
 
 #Collect contigs
 
 rpc <- rps$contig
 rdc <- rpd$contig
 mbc <- mbs$contig
 
 
 cnts <- unique(c(rpc, rdc, mbc))
 
 #Extract plasmidic contigs

 nms <- names(mtg)
 idx <- which(nms %in% cnts)
 
 hit <- mtg[idx]
 
 
 #Parsing output tables
 
 names(rps) <- c("Rep_ORF", "Rep_len", "Rep_type", "Rep_score", "contig")
 names(mbs) <- c("Mob_ORF", "Mob_len", "MOB_group", "MOB_score", "alifrom", "alito", "contig")
 names(rpd) <- c("Rep_type", "contig", "Rep_ORF")
 
 rpd1 <- rpd %>%
         group_by(contig) %>%
         summarise_all(funs(paste(., collapse = ',')))     
    
 rps1 <- rps %>% 
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 mbs1 <- mbs %>%
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 ftb <- full_join(rps1, mbs1, by = "contig")
 
 ftb1 <- full_join(ftb, rpd1, by = "contig")
 
 fct <- ftb1$contig
 nht <- names(hit)
 len <- numeric(0)
 
 for (i in 1:length(fct)){
   
   ct <- fct[i]
   idx <- which(nht == ct)
   cnl <- width(hit[idx])
   
   len <- c(len, cnl)
   
 }
 
 ftb1$contig_length <- len
 
 ftb2 <- ftb1[,c(c("contig","Rep_type.x","MOB_group","Rep_type.y","contig_length"))]
 colnames(ftb2)<- c("Contig", "Inc_group", "MOB_group", "RIP_domain", "contig_length")
 
 #Writing final results
 
 writeXStringSet(hit, "Plasmids_contigs.fasta")
 write_delim(ftb2, "Plasmid_Report.tsv", delim = "\t")
 
 