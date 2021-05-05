#!/bin/bash/env Rscript

#Rscript Retrieve_RIP_plasmids.R Filtered_classif.tsv Rep_domains.tsv Mob_table.tsv assembly.fa

 args = commandArgs(trailingOnly=TRUE)

 tbr = args[1]
 tbd = args[2]
 tbm = args[3]
 asm = args[4]

 tbr = 'Filtered_classif.tsv'
 tbd = 'Rep_domains.tsv'
 tbm = 'Mob_table.tsv'
 asm = 'assembly.fa'
 
 library(tidyverse)
 library(Biostrings)
 
 #Load data
 
 rps <- read_delim(tbr, delim="\t")
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
 
 if ( nrow(rpd) > 0 ) {
 
 names(rpd) <- c("Rep_domain", "contig", "Rep_domain_ORF")
 rpd1 <- rpd %>%
         group_by(contig) %>%
         summarise_all(funs(paste(., collapse = ',')))     
 } 
 
 if ( nrow(rps) > 0 ) {
 
 names(rps) <- c("Rep_type", "Rep_ORF", "Rep_score", "Rep_length", "contig")   
 rps1 <- rps %>% 
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 }
 
 if ( nrow(mbs) > 0 ) {
 
 names(mbs) <- c("Mob_ORF", "Mob_len", "MOB_group", "MOB_score", "alifrom", "alito", "contig")
 mbs1 <- mbs %>%
         group_by(contig) %>% 
         summarise_all(funs(paste(., collapse = ',')))
 
 }
 
 
 
 
 if ( length(mbc)>0 & length(rpc)>0 & length(rdc)>0 ) {
     
 mtb  <- full_join(rpd1, rps1, by = "contig")
 mtb1 <- full_join(mtb, mbs1, by = "contig")
     
 ftb2 <- mtb1[,c(c("contig","Rep_domain","MOB_group","Rep_type"))]
 colnames(ftb2)<- c("Contig", "RIP_domain", "MOB_group", "Inc_group")
 
 } else if ( length(mbc)== 0 & length(rpc)>0 & length(rdc)>0 ) {
     
 mtb1  <- full_join(rpd1, rps1, by = "contig")
     
 ftb2 <- mtb1[,c(c("contig","Rep_domain","Rep_type"))] 
 colnames(ftb2)<-  c("Contig", "RIP_domain", "Inc_group")

 ftb2$MOB_group <- rep(NA, length(ftb2$Contig))
 
 } else if ( length(mbc) > 0 & length(rpc)==0 & length(rdc)>0 ) {
 
 mtb1  <- full_join(rpd1, mbs1, by = "contig")
     
 ftb2 <- mtb1[,c(c("contig","Rep_domain","MOB_group"))]
 colnames(ftb2)<-  c("Contig", "RIP_domain", "MOB_group")
 ftb2$Inc_group <- rep(NA, length(ftb2$Contig))
 
 } else if ( length(mbc) > 0 & length(rpc)>0 & length(rdc)==0 ) {
    
 mtb1  <- full_join(rps1, mbs1, by = "contig")
     
 ftb2 <- mtb1[,c(c("contig","MOB_group","Rep_type"))]
 colnames(ftb2)<-  c("Contig", "MOB_group", "Inc_group")
 ftb2$RIP_domain <- rep(NA, length(ftb2$Contig))
 
 } else if ( length(mbc)==0 & length(rpc)==0 & length(rdc)>0 ) {
 
 ftb2 <- rpd1[,c(c("contig","Rep_domain"))]
 colnames(ftb2)<-  c("Contig", "RIP_domain")
 ftb2$Inc_group <- rep(NA, length(ftb2$Contig))
 ftb2$MOB_group <- rep(NA, length(ftb2$Contig))
    
 } else if ( length(mbc)==0 & length(rpc)>0 & length(rdc)==0 ) {
    
 ftb2 <- rps1[,c(c("contig","Rep_type"))]
 colnames(ftb2)<-  c("Contig", "Inc_group")
 ftb2$RIP_domain <- rep(NA, length(ftb2$Contig))
 ftb2$MOB_group <- rep(NA, length(ftb2$Contig))
 
 } else if ( length(mbc)>0 & length(rpc)==0 & length(rdc)==0 ) {
    
 ftb2 <- mbs1[,c(c("contig","MOB_group"))]
 colnames(ftb2)<-  c("Contig", "MOB_group")
 ftb2$RIP_domain <- rep(NA, length(ftb2$Contig))
 ftb2$Inc_group <- rep(NA, length(ftb2$Contig))
    
 }
 
 
 fct <- ftb2$Contig
 nht <- names(hit)
 len <- numeric(0)
 
 for (i in 1:length(fct)){
     
     ct <- fct[i]
     idx <- which(nht == ct)
     cnl <- width(hit[idx])
     
     len <- c(len, cnl)
     
 }
 
 ftb2$contig_length <- len
 
 #Writing final results
 
 writeXStringSet(hit, "Plasmids_contigs.fasta")
 write_delim(ftb2, "Plasmid_Report.tsv", delim = "\t")
 
 
