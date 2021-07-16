#!/usr/bin/env Rscript

 args = commandArgs(trailingOnly=TRUE)

 teb = args[1]
 
 
 library(tidyverse)
 

 tab <- read.table(teb, 
                  header = FALSE,
                  sep = "", 
                  blank.lines.skip = TRUE,
                  skipNul = TRUE,
                  col.names= c("RIP","taccession","tlen","queryname", "qaccession","qlen","Evalue",
                              "score","bias","num","of","cEvalue","iEvalue", "domscore","dombias",
                              "hmmfrom","hmmto","alifrom","alito","from_env","to_env","env"))

  tib <- as_tibble(tab)

  tsd <- tibble('RIP'        = character(),
                'tlen'       = character(),
                'queryname'  = character(),
                'qaccession' = character(),
                'qlen'       = numeric(),
                'Evalue'     = numeric(),
                'score'      = numeric(),
                'hmmfrom'    = numeric(),
                'hmmto'      = numeric(),
                'alifrom'    = numeric(),
                'alito'      = numeric())
 
  RIP <- unique(tib$RIP)
  
  #Name of RIPs in list l1 
  nl1 <- character(0)
  l1 <- list()
  x <- 0
 
  for (i in 1:length(RIP)) {
    
    ri <- as.character(RIP[i])
    
    idx <- which(tib$RIP == ri)
    tab1 <- tib[idx,]
    nls <- nrow(tab1)
    
    
    if (nls > 1) {
      
      nls <- c(nls, RIP)
      
      x <- x + 1
      
      
      e1 <- tab1[1,]$alito
      s2 <- tab1[2,]$alifrom
      
      if (is.na(s2)==TRUE){
        dom1a <- as.character(tab1[1,4])
        dom1  <- as.vector(c(dom1a,"no_2nd_domain"))
      } else {
        if (e1<s2){
          dom1a <- as.character(tab1[1,4])
          dom1  <- as.vector(c(dom1a, "not_over"))
        } else {
          bs1 <- as.numeric(tab1[1,]$score)
          bs2 <- as.numeric(tab1[2,]$score)
          if (bs1<bs2){
            dom1a <- as.character(tab1[2,]$queryname)
            dom1  <- as.vector(c(dom1a,"overlapped"))
          } else {
            dom1a <- as.character(tab1[1,]$queryname)
            dom1  <- as.vector(c(dom1a,"overlapped"))
          }
        }
      }
      
      
      
      e2 <- tab1[2,]$alito
      s3 <- tab1[3,]$alifrom
      
      if (is.na(s3)==TRUE){
        dom2a <- as.character(tab1[2,]$queryname)
        dom2  <- as.vector(c(dom2a,"no_3rd_domain"))
      } else {
        if (e2<s3){
          dom2a <- as.character(tab1[2,]$queryname)
          dom2  <- as.vector(c(dom2a, "not_over"))
        } else {
          bs2 <- as.numeric(tab1[2,]$score)
          bs3 <- as.numeric(tab1[3,]$score)
          if (bs2<bs3){
            dom2a <- as.character(tab1[3,]$queryname)
            dom2  <- as.vector(c(dom2a,"overlapped"))
          } else {
            dom2a <- as.character(tab1[2,]$queryname)
            dom2  <- as.vector(c(dom2a,"overlapped"))
          }
        }
      }
      
      e3 <- tab1[3,]$alito
      s4 <- tab1[4,]$alifrom
      
      if (is.na(s4)==TRUE){
        dom3a <- as.character(tab1[3,]$queryname)
        dom3  <- as.vector(c(dom3a,"no_4th_domain"))
      } else {
        if (e3<s4){
          dom3a <- as.character(tab1[3,]$queryname)
          dom3  <- as.vector(c(dom3a, "not_over"))
        } else {
          bs3 <- as.numeric(tab1[3,]$score)
          bs4 <- as.numeric(tab1[4,]$score)
          if (bs3<bs4){
            dom3a <- as.character(tab1[4,]$queryname)
            dom3  <- as.vector(c(dom3a,"overlapped"))
          } else {
            dom3a <- as.character(tab1[3,]$queryname)
            dom3  <- as.vector(c(dom3a,"overlapped"))
          }
        }
      }
      
      
      e4 <- tab1[4,]$alito
      s5 <- tab1[5,]$alifrom
      
      if (is.na(s5)==TRUE){
        dom4a <- as.character(tab1[4,]$queryname)
        dom4  <- as.vector(c(dom4a,"no_5th_domain"))
      } else {
        if (e4<s5){
          dom4a <- as.character(tab1[4,]$queryname)
          dom4  <- as.vector(c(dom4a, "not_over"))
        } else {
          bs4 <- as.numeric(tab1[4,]$score)
          bs5 <- as.numeric(tab1[5,]$score)
          if (bs4<bs5){
            dom4a <- as.character(tab1[5,]$queryname)
            dom4  <- as.vector(c(dom4a,"overlapped"))
          } else {
            dom4a <- as.character(tab1[4,]$queryname)
            dom4  <- as.vector(c(dom4a,"overlapped"))
          }
        }
      }
      
      
      Arq <- as.vector(c(dom1,dom2,dom3,dom4))  
      Arq <- Arq[!is.na(Arq)]
      
      l1[[x]] <- Arq
      nl1 <- c(nl1, as.character(ri))
      
      
    } else {
      
      slc <-  c(1, 3, 4, 5, 6, 7, 8, 16, 17, 18, 19)
      
      stb <- tab1[1,slc]
      
      tsd <- rbind(tsd, stb)
      
      
    }
    
  }
  
  
  saveRDS(l1, "Domain_Architecture.RDS")
  nl2 <- as.tibble(nl1)
  write_delim(nl2, "multi_dom_RIP.tsv")
  write_delim(tsd, "single_dom_RIP.tsv")
  
  
