#!/bin/bash/env Rscript

args = commandArgs(trailingOnly=TRUE)

paf = args[1]         # Input paf table
fst = args[2]         # Input fasta file
rds = args[3]         # Plsdb length and names

library(tidyverse)
library(Biostrings)
library(readr)

tab.all <- read.csv(paf,sep='\t',header=F)

colnames(tab.all) <- c('qid','qlen','qst','qend','strand','sid','slen','sst','send','match','len','qual')

fs1   <- readDNAStringSet(fst)
nouts <- unique(names(fs1))

cns <- strsplit(nouts[1], split = "_")[[1]][1]

ctg <- character(0)
sbj <- character(0)
similarities <- list()

sbt <- character(0)

for (x in 1:length(nouts)) {

  # Parse PAF table
  tnc  <- nouts[x]
  tab  <- tab.all[which(tab.all$qid == tnc),]
  tab1 <- tab[1,]

  pid  <- (tab1$match/tab1$len)*100
  qcov <- (tab1$len/tab1$qlen)*100
  S    <- (pid*qcov)/100

   if (dim(tab)[1]==0) {

    similarities[[x]] <- 0

  } else {

     similarities[[x]] <- S

    dtf <- as.data.frame(table(tab$sid))
    csx <- dtf[order(dtf$Freq, decreasing = TRUE),]
    dbh <- csx$Var1[1]

    tq  <- tab$qid
    dbq <- tq[1]

    sbj[x]	<- as.character(as.vector(dbh))
  }

    sbt <- c(sbt, as.character(tab1$sid))

}


# Filtering mapped contigs id.qcov > 45 for plasmids

mean.vect <- rapply(similarities, mean)

sim       <- mean(mean.vect)

#collecting data of matching plasmid

plsdb_table <- readRDS(rds)
plsdb_names <- plsdb_table$plasmid_name
plsdb_len   <- plsdb_table$plasmid_length

sbl <- as.data.frame(table(sbt))
plf <- sbl[order(sbl$Freq, decreasing = TRUE),]
pl1 <- as.character(plf$sbt[1])
spl <- grep(pl1, plsdb_names)

psn <- plsdb_names[spl]
pll <- plsdb_len[spl]

# Output data

df <- tibble('Contig_name' = character(),'S-distance' = numeric(), 'Ref_name' = character(), 'Ref_length' = numeric())

dat <- c(cns, sim, psn, pll)


 if (sim >= 45) {
     
     df[1,1]   <- cns
     df[1,2]   <- sim
     df[1,3]   <- psn
     df[1,4]   <- pll
     
     write_delim(df, delim = "\t", path = paste0(cns,".tsv"), col_names = FALSE)

 } else {
   
      
      df[1,1]   <- cns
      df[1,2]   <- sim
      df[1,3]   <- NA
      df[1,4]   <- NA
   
      write_delim(df, delim = "\t", path = paste0(cns,".tsv"), col_names = FALSE) 
   
 }


