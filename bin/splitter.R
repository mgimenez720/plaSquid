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

q <- readDNAStringSet(query)
l <- length(q)

qnames       <- c()
lst <- list(0)


for (x in 1:l) {
  
  
  s  <- as.list(as.character(q[x]))
  n  <- names(s)
  n2 <- strsplit(n,' ')[[1]][1]
  v  <- lapply(s,s2c)[[1]]
  
  qnames <- c(qnames,n2)
  
  spl <- spliter(w=fsize,s=wsize,v=v)
  len <- length(spl)
  nms <- paste(n2,1:len,sep='_')
  
  names(spl) <- nms
  
  fout <- tempfile(tmpdir='.',fileext='.fasta')
  
  write.fasta(spl, names = names(spl), file.out=fout)
  
  #cmnd <- paste("cat", fout, ">> plasmid.split", sep = " ")
  #system(cmnd)
  
  #rmr <- paste("rm -r ", fout, sep="")
  #system(rmr, ignore.stdout = TRUE)
}


