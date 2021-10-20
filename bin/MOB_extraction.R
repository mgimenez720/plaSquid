#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

 library(Biostrings)
 library(tidyverse)

 prt = args[1]
 tbd = args[2]
 
 
 prots <- readAAStringSet(prt)
 npr <- names(prots)
 
 nms <- character(0)
 
 for(i in 1:length(npr)) {
   
   nri <- npr[i]
   
   nm <- strsplit(nri, split = " ")[[1]][1]
   
   nms <- c(nms, nm)
 }
 
 
 tab   <- read_delim(tbd, delim = "\t")
 tmd   <- tab$Mob_det
 idx   <- which(nms %in% tmd)
 
 seqs <- prots[idx]
 
 nmf <- character(0) 
 
 for(i in 1:length(tab$Mob_det)) {
   
   tbi <- tab$Mob_det[i]
   mbi <- tab$query_name[i]
   nmi <- paste0(tbi,"_",mbi)
   
   j <- grep(tbi,names(seqs))
   names(seqs)[j] <- nmi 
   
 } 
 
 
 writeXStringSet(seqs, "MOB_seqs.faa")
 
 