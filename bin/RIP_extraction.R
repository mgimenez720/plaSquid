#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

prt = args[1] #"prots.faa"
inc = args[2] #"Filtered_Classif.tsv"
dom = args[3] #"Rep_domains.tsv"


library("Biostrings")
library("tidyverse")

orf <- readAAStringSet(prt)
nrf <- names(orf)

nms <- character(0)

for(i in 1:length(nrf)) {
  
  nri <- nrf[i]
  
  nm <- strsplit(nri, split = " ")[[1]][1]
  
  nms <- c(nms, nm)
}

tbi  <- read_delim(inc, delim = "\t")
idd  <- grep("_", tbi$query_name)
rpi  <- tbi$query_name[idd]
tbi1 <- tbi[idd,]


tbd <- read_delim(dom, delim = "\t")
idx <- which(!is.na(tbd$Rep_ORF))
rpo <- tbd$Rep_ORF[idx]

tbi2 <- tbi1[,c("contig","query_name","Inc_det")]
tbd1 <- tbd[,c("contig","Rep_ORF","Rep_type")]

colnames(tbi2) <- colnames(tbd1)

tbs <- rbind(tbd1, tbi2)

tbb <- tbs[which(!is.na(tbs$contig)),]

tbf <-tbb[order(tbb$contig, tbb$Rep_ORF),]

rip <- unique(tbf$Rep_ORF)
idf <- which(nms %in% rip)
rep <- orf[idf]


tbf1 <- tbf %>% group_by_at(vars(Rep_ORF)) %>%
         summarize_all(paste, collapse="#")

nmf <- character(0)

for(i in 1:length(tbf1$Rep_ORF)){
  
  rpi <- tbf1$Rep_ORF[i]
  rti <- tbf1$Rep_type[i]
  
  nmi <- paste(rpi,rti, sep ="#")
  
  nmf <- c(nmf, nmi)
  
}

names(rep) <- nmf

writeXStringSet(rep, "RIP_seqs.faa")
