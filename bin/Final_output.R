#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

frs = args[1] #Result.tsv
dna = args[2] #Result.fasta
bct = args[3] #Bact_120.tsv


library(tidyverse)
library(Biostrings)

   

   dn1 <- readDNAStringSet(dna)

     rst <- as_tibble(read.table(frs, header=TRUE, sep = "\t", 
                                 comment.char = '#', 
                                 col.names= c("Contig", "name", "Sim_dist", "plsdb_match","Match_length",
                                              "RIP_domain", "MOB_group", "Rep_type", "Contig_length"),
                                 fill = TRUE, stringsAsFactors = FALSE))
     
     names(dn1) <- rst$Contig

     hmm <- as_tibble(read.table(bct, 
                                 header=FALSE, 
                                 sep = "", 
                                 comment.char = '#', 
                                 col.names= c("query_name","taccession","tlen","Marker_gene",
                                              "qaccession","qlen","Evalue","score","bias",
                                              "num","of","CEvalue","iEvalue","domscore","dombias",
                                              "hmmfrom","hmmto","alifrom","alito","envfrom",
                                              "envto","acc"), 
                                 fill = TRUE, 
                                 stringsAsFactors = FALSE))

     
     if (nrow(hmm) >= 1) {
     
      cnt <- character(0)

         for (i in 1:length(hmm$query_name)){
  
                    ri  <- hmm$query_name[i]
                    rc  <- strsplit(ri, split = "_")[[1]][1]
                    cnt <- c(cnt, rc)
            }
  
      ucn <- unique(cnt)

   
   
     idf <-  which(!rst$Contig %in% ucn)
      
      
     fdn <- dn1[idf]   
     frt <- rst[idf,]
     
     names(fdn) <- frt$name
          
     writeXStringSet(fdn, "plaSquid_result.fasta")
     write_delim(frt, "plaSquid_result.tsv", delim = "\t")
     
     } else {

       names(dn1) <- rst$name       
       
       writeXStringSet(dn1, "plaSquid_result.fasta")
       write_delim(rst, "plaSquid_result.tsv", delim = "\t")   
       
       
     }
     
