#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

query = args[1]         # Input file
fsize = 1500            # Size of fragments to be generated from contigs 
wsize = 250             # Sliding window size


library(Biostrings,quietly=T)
library(seqinr,quietly=T)

spliter <- function(
  
  w=win,
  s=winst,
  v
  
) {
  
  n  <- length(v)
  x  <- 1
  x2 <- w
  o  <- 1
  l  <- list()
  
  if ( n>w & x2<=(n-s) ) {
    
    while ( x2<=(n-s) ) {
      
      ini    <- x
      fin    <- x2
      l[[o]] <- v[ini:fin]
      x      <- x+s
      x2     <- x2+s
      o      <- o+1
      
    }
    
    l[[o]] <- v[(fin+1):n]
    
  } else {
    
    l[[1]] <- v
  }
  
  return(l)
}     

#------#

#Rename contigs and move files

fst  <- readDNAStringSet(query)
ctr  <- 0
nm.s <- character(0)

for (x in 1:length(fst)) {
  
  fx <- fst[x]
  
  s  <- as.list(as.character(fst[x]))
  v  <- lapply(s,s2c)[[1]]
  
  ctr  <- ctr+1
  nm   <- paste0("Contig-",ctr)
  nmf  <- paste0(nm,".fst")
  nm.s <- c(nm.s, nm)
  
  spl <- spliter(w=fsize,s=wsize,v=v)
  len <- length(spl)
  nms <- paste(nm,1:len,sep='_')
  
  names(spl) <- nms
  
  write.fasta(spl, names = names(spl), file.out=nmf)

}



