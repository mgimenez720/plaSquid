#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

hmt = args[1]


library(tidyverse)


tab <- as_tibble(read.table(hmt, header=FALSE, sep = "", comment.char = '#', 
                            col.names= c("Mob_det","taccession","tlen","query_name",
                                         "qaccession","qlen","Evalue","score","bias",
                                         "num","of","CEvalue","iEvalue","domscore","dombias",
                                         "hmmfrom","hmmto","alifrom","alito","envfrom",
                                         "envto","acc", "DoT"), fill = TRUE, stringsAsFactors = FALSE))


mbs <- c("MOBB", "MOBC", "MOBF", "MOBH", "MOBP1", "MOBP2", "MOBP3", "MOBQ", "MOBT", "MOBV", "MOBM")

smb <- c(92.7, 96.6, 332.2, 81.1, 74.5, 480, 105, 61, 115.9, 68, 584)

tbf <- tibble(Mob_det    = character(0),
              tlen       = integer(length = 0L),
              query_name = character(0),
              score      = numeric(),
              alifrom    = integer(length = 0L),
              alito      = integer(length = 0L))


for (i in 1:length(mbs)) {
  
  mb <- mbs[i]
  sc <- smb[i]
  
  tbi <- subset.data.frame(tab, tab$query_name == mb, select = c(Mob_det, tlen, query_name, score, alifrom, alito))
  
  tbh <- subset.data.frame(tbi, tbi$score >= sc)
  
  tbu <- distinct(tbh, Mob_det, query_name, .keep_all = TRUE)
  
  tbf <- rbind(tbf, tbu)
  
}
 
if (nrow(tbf)>0) {

 cn <- character(0)

 for (i in 1:length(tbf$Mob_det)){
  
  pri <- tbf$Mob_det[i]
  cnt <- strsplit(pri, split = "_")[[1]][1]
   
  cn <- c(cn, cnt) 
  
 }
 
 tbf$contig <- cn

}
 write_delim(tbf, "Mob_table.tsv", delim = "\t")

