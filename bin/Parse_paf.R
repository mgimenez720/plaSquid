#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

paf = args[1]         # Input paf table
fst = args[2]         # Input fasta file
cnt = args[3]         # Splitted fasta file 
rds = args[4]         # Plsdb length and names

library(tidyverse)
library(Biostrings)
library(readr)

if(file.size(paf) > 0) {
  
tab.all <- read.csv(paf,sep='\t',header=F)

colnames(tab.all) <- c('qid','qlen','qst','qend','strand','sid','slen','sst','send','match','len','qual')

fs1   <- readDNAStringSet(fst)
nouts <- names(fs1)
louts <- width(fs1)
nct <- character(0)

tig <- readDNAStringSet(cnt)
ntg <- names(tig)
rm(tig)

for (i in 1:length(fs1)) {
  
  nout <- nouts[i]
  cns  <- strsplit(nout, split = " ")[[1]][1]
  cnx  <- paste0(cns,"_")
  nct  <- c(nct, cnx)
  
}

ctg <- character(0)
sbj <- character(0)
similarities <- list()

sbt <- character(0)

for (x in 1:length(nct)) {

  # Parse PAF table
  tnc  <- nct[x]
  idx  <- grep(tnc, tab.all$qid, fixed = TRUE)
  tab  <- tab.all[idx,]
  
  sbs <- unique(grep(tnc, ntg, value = TRUE))
  
  
  if (nrow(tab) > 0) {
  
  S <- numeric(length = 0L)
  
  for (i in 1:length(sbs)) {
    
  sb <- sbs[i]  
  
  tab2 <- subset.data.frame(tab, tab$qid == sb)
  
  if (nrow(tab2)>0) {
  
  tab1 <- tab2[1,]

  pid  <- (tab1$match/tab1$len)*100
  qcov <- (tab1$len/tab1$qlen)*100
  s1   <- (pid*qcov)/100

  S <- c(S, s1)
  
  } else {
    
  S <- c(S, 0) 
    
  }
  }
  
  similarities[[x]] <- S
  
  dtf <- as.data.frame(table(tab$sid))
  csx <- dtf[order(dtf$Freq, decreasing = TRUE),]
  dbh <- csx$Var1[1]
  
  tq  <- tab$qid
  dbq <- tq[1]
  
  sbj[x]	<- as.character(as.vector(dbh))
  
    
  } else {

    similarities[[x]] <- 0

    sbj[x] <- NA
  
  }

   
  
}

sim <- rapply(similarities, mean)

#collecting data of matching plasmid

plsdb_table <- readRDS(rds)
plsdb_names <- plsdb_table$plasmid_name
plsdb_len   <- plsdb_table$plasmid_length

# Filtering mapped contigs id.qcov > 45 for plasmids

psn <- character(0)
psl <- numeric(length=0L)
idx <- numeric(length=0L)

for (j in 1:length(nouts)){

if (sim[j] > 50) {  
    
  sbi <- sbj[j]
  spl <- grep(sbi, plsdb_names, fixed = TRUE)
  ps1 <- plsdb_names[spl]
  pl1 <- plsdb_len[spl]
  idx <- c(idx, j)
  
} else {
  
  ps1 <- NA
  pl1 <- NA
  
}
  
  psn <- c(psn, ps1)
  psl <- c(psl, pl1)

}


# Output data

  df <- tibble('Contig_name'   = nouts,
               'S-distance'    = sim,
               'Ref_name'      = psn, 
               'Ref_length'    = psl,
               'Contig_length' = louts)


} else {
  
  df <- tibble('Contig_name'   = NA,
               'S-distance'    = NA,
               'Ref_name'      = NA, 
               'Ref_length'    = NA,
               'Contig_length' = NA)
  
  
} 

     df1 <- subset.data.frame(df, df$`S-distance`>50)
     
     fsj <- fs1[idx]
     writeXStringSet(fsj, "Minidist_plasmid.fasta")
     write_delim(df1, delim = "\t", path = "Minidist_result.tsv", col_names = TRUE)

 

