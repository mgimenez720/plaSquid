#!/usr/bin/Rscript


 args = commandArgs(trailingOnly=TRUE)

 trb = args[1]
 
 library(tidyverse)

 
 tab <- as_tibble(read.table(trb, header = TRUE, stringsAsFactors = FALSE))
   
   
# Filter IncFII and IncZ plasmids
   
   pls <- unique(tab$contig)
   nd  <- character(0)
   
   for (i in 1:length(pls)){
     
     pli <- pls[i]
     
     tbp <- subset.data.frame(tab, tab$contig==pli)
     
     idx <- grep("IncFII", tbp$query_name, fixed = TRUE)
     
     if (length(idx) == 0) {
       
       next
       
     } else {
       
       cnt  <- tbp$contig
       
       tbp1 <- subset.data.frame(tbp, tbp$contig==cnt)
       
       qnt  <- tbp1$query_name
       
       xdi  <- grep("IncZ", qnt, fixed = TRUE)
       
       if (length(xdi) > 0) {
         
         prot <- as.vector(tbp1$Inc_det[xdi])
         
         nd   <- c(nd, prot)
         
       } else {
         
         next
         
       }
     }
   }
   
   tbi <- tab$Inc_det
   xdd <- which(tbi %in% nd)

   if (length(xdd) > 0) {
      
    atb <- tab[-xdd,] 
   
     } else {
    
    atb <- tab       
      
   }
   
# Filter Col and Gram-positive plasmids
   
   pls <- unique(atb$contig)
   nd  <- character(0)
   
   for (i in 1:length(pls)){
     
     pli <- pls[i]
     
     tbp <- subset.data.frame(atb, atb$contig==pli)
     
     idx <- grep("Col", tbp$query_name, fixed = TRUE)
     
     if (length(idx) == 0) {
       
       next
       
     } else {
       
       cnt  <- tbp$contig
       
       tbp1 <- subset.data.frame(tbp, tbp$contig==cnt)
       
       qnt  <- tbp1$query_name
       
       xdi  <- grep("Inc13", qnt, fixed = TRUE)
       
       if (length(xdi) > 0) {
         
         prot <- as.vector(tbp1$Inc_det[xdi])
         
         nd   <- c(nd, prot)
         
       } else {
         
         next
         
       }
     }
   }
   
   
   tbi <- atb$Inc_det
   xdd <- which(tbi %in% nd)
   
   if (length(xdd) > 0) {
      
      tib <- atb[-xdd,] 
      
   } else {
      
      tib <- atb       
      
   }
   
# Write final table
   
   write_delim(tib, "Filtered_Classif.tsv", delim = "\t")
   