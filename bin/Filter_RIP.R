#!/bin/bash/env Rscript

    args = commandArgs(trailingOnly=TRUE)
    
    mdt = args[1]  #multi_domain_tab
    sdt = args[2]  #single_domain_tab
    rds = args[3]  #RIP_candidates_Arch_list
    arq = args[4]  #RIP_Arch_list
    
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
    dmn <- character(0)
    
    for (x in 1:length(lsds)) {
      
      lsd   <- lsds[x]
  
      s_dom <- subset.data.frame(tbs, tbs$queryname == lsd)
      s_dom$score <- as.numeric(as.character(s_dom$score))
      
      if ( nrow(s_dom) == 0) {
        
        next
      
        } else {
        
        
        tsh   <- s_dom[(s_dom$score > as.numeric(lsdvs[lsd])) ,]
        hit   <- as.vector(tsh$RIP)
        dum   <- rep(lsd, length(hit))
        
        sd1 <- c(sd1,hit)
        dmn <- c(dmn, dum)
      }
      
    }
    
    c1t <- character(0)
    
    if (length(sd1) > 0) {
    
    for(i in 1:length(sd1)){
      
      ht <- sd1[i]
      cn <- strsplit(ht, split = "_")[[1]][1]
      c1t <- c(c1t, cn) 
    
    }
    
    ssp <- tibble("Rep_type" = dmn,
                  "contig"   = c1t,
                  "Rep_ORF"  = sd1)
    
    } else {
      
    ssp <- tibble("Rep_type" = NA,
                  "contig"   = NA,
                  "Rep_ORF"  = NA) 
      
    }
    
    # Filtering single-domain RIPs by bit-score and length.
    
    Sdm <- c("PriCT_1","Rep_1","Rep_2","Rep_3","RepL","Rep_trans")
    
    sdom        <- subset.data.frame(tbs, tbs$queryname == "PriCT_1")
    PriCT_1     <- sdom[(sdom$score > 49 & sdom$tlen < 500 & sdom$tlen > 420) ,]
    PriCT_1_sdl <- as.vector(PriCT_1$RIP)
    rpr         <- rep("PriCT_1", length(PriCT_1_sdl))
    
    sdom      <- subset.data.frame(tbs, tbs$queryname ==  "Rep_1")
    Rep_1A    <- sdom[(sdom$score > 37 & sdom$tlen > 130) ,]
    Rep_1B    <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_1     <- rbind(Rep_1A, Rep_1B)
    Rep_1_sdl <- as.vector(Rep_1$RIP)
    rp1       <- rep("Rep_1", length(Rep_1_sdl))
    
    
    sdom      <- subset.data.frame(tbs, tbs$queryname ==  "Rep_3")
    rep3      <- sdom[(sdom$score > 45),]
    Rep_3_sdl <- as.vector(rep3$RIP)
    rp3       <- rep("Rep_3", length(Rep_3_sdl))
    
    sdom      <- subset.data.frame(tbs, tbs$queryname == "RepL")
    RepL      <- sdom[(sdom$score > 85 & sdom$tlen > 90) ,]
    RepL_sdl  <- as.vector(RepL$RIP)
    rpl       <- rep("RepL", length(RepL_sdl))
    
    sdom          <- subset.data.frame(tbs, tbs$queryname == "Rep_trans")
    Rep_trans     <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_trans_sdl <- as.vector(Rep_trans$RIP)
    rpt           <- rep("Rep_trans", length(Rep_trans_sdl))
    
    
    hts <- c(PriCT_1_sdl,  Rep_1_sdl, Rep_3_sdl, RepL_sdl, Rep_trans_sdl)
    sdo <- c(rpr, rp1, rp3, rpl, rpt)
    
    gtc <- character(0)
    
    if (length(hts) > 0 ) {
    
    for (i in 1:length(hts)){
      
      tih <- hts[i]
      itc <- strsplit(tih, split = "_")[[1]][1]
      
      gtc <- c(gtc, itc)
    
    }
    
    sop <- tibble("Rep_type"  = sdo,
                  "contig"    = gtc,
                  "Rep_ORF"   = hts)
    
    } else {
      
    sop <- tibble("Rep_type"  = NA,
                  "contig"    = NA,
                  "Rep_ORF"   = NA)
      
    }
    
    # Multi-domain RIP candidate filtering
   
    idx <- which(lsr %in% ar1)
    htm <- as.character(tbm[idx,])
    
    cnn <- character(0)
    
    if (length(htm > 0)) {
    
    for (i in 1:length(htm)){
      
      hti <- htm[i]
      ctg <- strsplit(hti, split = "_")[[1]][1]
      cnn <- c(cnn, ctg)
      
    }
    
      nmm <- length(htm)
      cda <- rep("Conserved Domain Arch", nmm)
    
      mop <- tibble("Rep_type" = cda,
                    "contig"  = cnn,
                    "Rep_ORF" = htm)
    
    } else {
      
      mop <- tibble("Rep_type" = NA,
                    "contig"   = NA,
                    "Rep_ORF"  = NA)
      
      
    }
    
      
   ftb <-  rbind(ssp, sop, mop)
    
   write_delim(ftb, "Rep_domains.tsv", delim = "\t")
   
   
   
