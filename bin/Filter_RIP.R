#!/bin/bash/env Rscript

    args = commandArgs(trailingOnly=TRUE)
    
    mdt = args[1]
    sdt = args[2]
    rds = args[3]
    arq = args[4]
    
    library(tidyverse)
    
    #Read data
    tbm <- read.table(mdt, header=TRUE, sep = "\t")
    tbs <- read.table(sdt, header=TRUE, sep = " ")
    lsr <- readRDS(rds)
    ar1 <- readRDS(arq)
    
    #Filter single domain RIP candidates
    
    #Filtering single-domain RIPs by bit-score.
    
    lsds  <- c("IncFII_repA","RepA_C", "RepA_N", "RepC", "Replicase", "Rop", 
               "RPA", "RP-C", "TrfA", "Bac_RepA_C", "RepB-RCR_reg", "RP-C_C")
    
    lsdvs <- c(10.0, 45, 38, 76, 76, 77.1, 67.8, 45, 87, 30, 24, 34)
    
    names(lsdvs) <- lsds
    
    sd1 <- character(0)
    
    for (x in 1:length(lsds)) {
      
      lsd   <- lsds[x]
  
      s_dom <- subset.data.frame(tbs, tbs$queryname == lsd)
      s_dom$score <- as.numeric(as.character(s_dom$score))
      
      if ( nrow(s_dom) == 0) {
        
        next
      
        } else {
        
        
        tsh   <- s_dom[(s_dom$score > as.numeric(lsdvs[lsd])) ,]
        hit   <- as.vector(tsh$RIP)
        
        sd1 <- c(sd1,hit)
      
      }
      
    }
    
    
    # Filtering single-domain RIPs by bit-score and length.
    
    Sdm <- c("PriCT_1","Rep_1","Rep_2","Rep_3","RepL","Rep_trans")
    
    sdom        <- subset.data.frame(tbs, tbs$queryname == "PriCT_1")
    PriCT_1     <- sdom[(sdom$score > 49 & sdom$tlen < 500 & sdom$tlen > 420) ,]
    PriCT_1_sdl <- as.vector(PriCT_1$RIP)
    
    sdom      <- subset.data.frame(tbs, tbs$queryname ==  "Rep_1")
    Rep_1A    <- sdom[(sdom$score > 37 & sdom$tlen > 130) ,]
    Rep_1B    <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_1     <- rbind(Rep_1A, Rep_1B)
    Rep_1_sdl <- as.vector(Rep_1$RIP)
    
    sdom      <- subset.data.frame(tbs, tbs$queryname == "Rep_2")
    Rep_2_sdl <- as.vector(sdom$RIP)
    
    sdom      <- subset.data.frame(tbs, tbs$queryname ==  "Rep_3")
    Rep_3_sdl <- as.vector(sdom$RIP)
    
    sdom     <- subset.data.frame(tbs, tbs$queryname == "RepL")
    RepL     <- sdom[(sdom$score > 85 & sdom$tlen > 90) ,]
    RepL_sdl <- as.vector(RepL$RIP)
    
    sdom          <- subset.data.frame(tbs, tbs$queryname == "Rep_trans")
    Rep_trans     <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_trans_sdl <- as.vector(Rep_trans$RIP)

    
    c()
    
    
    # Multi-domain RIP candidate filtering 
    
    cmp <- sapply(lst,identical,ar1)
    cmv <- as.vector(cmp)
    
    idx <- which(cmv == TRUE)
    
    htm <- tbm[idx,]
    
    