#!/usr/bin/Rscript


args = commandArgs(trailingOnly=TRUE)

hmt = args[1]
cmt = args[2]

library(tidyverse)


hmt <- as_tibble(read.table(hmt, header=FALSE, sep = "", comment.char = '#', 
            col.names= c("Inc_det","taccession","tlen","query_name",
                        "qaccession","qlen","Evalue","score","bias",
                        "num","of","CEvalue","iEvalue","domscore","dombias",
                        "hmmfrom","hmmto","alifrom","alito","envfrom",
                        "envto","acc"), fill = TRUE, stringsAsFactors = FALSE))

Incs <- c("Inc11_B", "Inc11", "Inc13_A", "Inc13_B", "Inc13_C", "Inc18", "Inc1", "Inc1_B",
          "Inc4_A", "Inc4B-9-10-14", "Inc7_A","Inc7_B","Inc8", "IncAC", "IncB-O-K-I", 
          "IncFI_RepB", "IncFI_RepE", "IncGU", "IncHIA", "IncHIB", "IncHI2", "IncLM",
          "IncN","IncP2","IncP7","IncP9","IncP","IncQ","IncR","IncT","IncW","IncX","IncZ")

sci <- c(350.2, 193.7, 491.0, 448.4, 530.2, 450.6, 532.6, 712.9, 568.4, 515.6, 489.0, 551.3,
         657.4, 655.7, 573.9, 438, 517.8, 860.4, 451.0, 567, 852.1, 466.2, 257.0, 392.1,
         674.8, 236.4, 593.5, 469.7, 658.9, 675.5, 449.4, 443.1, 515.0)

hdl <- tibble("Inc_det" = character(0),
              "tlen" = character(0),
              "query_name" = character(0),
              "score" = numeric(0))
       
 for (i in 1:length(Incs)) {
   
   inc <- Incs[i]
   scc <- sci[i]
   
   hm1 <- subset.data.frame(hmt, hmt$query_name == inc)
   hms <- subset.data.frame(hm1, subset = hm1$score >= scc, select = c("Inc_det", "tlen", "query_name", "score"))
   
   hdl <- rbind(hdl, hms)
   
 }
 
  cnt <- character(0)
 
  for (i in 1:length(hdl$Inc_det)){
    
    ri <- hdl$Inc_det[i]
    
    rc <- strsplit(ri, split = "_")[[1]][1]
    
    cnt <- c(cnt, rc)
  }

  hdl$contig <- cnt

# Threshold specific filtering of RNA incompatibility determinants

cm1 <- as_tibble(read.table(cmt, header=FALSE, sep = "", comment.char = '#', 
                           col.names= c("Inc_det","accession", "query_name","accession","mdl",
                                        "mdl_from","mdl_to","seq_from","seq_to","strand","trunc",
                                        "pass","gc","bias","score","E-value", "inc", "DoT"),
                           fill = TRUE, stringsAsFactors = FALSE))
         
 
 cmm <- c("Col156","Col3M","Col440I","Col440II","Col8282","Col(BS512)","ColE10","Col(IMGS31)","Col(IRGK)",
          "ColKP3","Col(MP18)","ColpVC","ColRNAI","Col(SD853)","Col(Ye4449)","Col(MG828)","IncFII_1_pKP91",
          "IncFII(29)_1_pUTI89","IncFII(p96A)_1_p96A","IncFII(pCoo)_1_pCoo","IncFII(pCRY)_1_pCRY","IncFII(pCTU2)_1_pCTU2",
          "IncFII(pECLA)_1_pECLA","IncFII(pHN7A8)_1_pHN7A8","IncFII(pKPX1)","IncFII(pMET)_1_pMET1","IncFII(pRSB107)_1_pRSB107",
          "IncFII(pSE11)_1_pSE11","IncFII(pYVa12790)_1_pYVa12790","IncFII(S)_1","IncFII(Y)_1_ps","IncFII_p14_Yersenia")


 scm <- c(138.2, 106.9, 115.4, 248.8, 193.6, 245.0, 179.6, 197.2,190.3, 298.4,182.9,178.0,113.2,173,188.8,182,
         210.5,266.6,587.6, 264,608.5,636.1,786.8,262.5,585.0,629.7,261.8,274.2,697.6,238.4,217.4,205.9)
 
 hdc <- tibble("Inc_det" = character(0),
               "query_name" = character(0),
               "score" = numeric(0))
 

 for (i in 1:length(cmm)) {
   
   
   cmi <- cmm[i]
   sci <- scm[i]
   
   cmk <- subset.data.frame(cm1, cm1$query_name == cmi)
   cmh <- subset.data.frame(cmk, cmk$score >= sci, select = c(Inc_det, query_name, score, mdl_from))
   
   if (cmi == "Col440I") {
     
     cmh1 <- subset.data.frame(cmh, cmh$mdl_from == 1, select = c(Inc_det, query_name, score))
     
   } else {
     
     cmh1 <- cmh[,1:3] 
     
   }
   
   hdc <- rbind(hdc, cmh1)
 }

  n <- length(hdc$Inc_det)
  hdc$tlen <- rep("NA", n)
  hdc$contig <- hdc$Inc_det
  

  ctb <- rbind(hdc, hdl)
  
  
  write_delim(ctb, "Classification_table.tsv", delim = "\t")

  
                           