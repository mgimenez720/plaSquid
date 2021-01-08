#!/usr/bin/env Rscript

#' \strong{minidist}: plasmidic contigs detection through database comparison
#'
#' @description this function looks for plasmid derived (meta)genomic contigs by mapping contigs subsequences
#' against a plasmidic database.  This function runs minimap2 for mapping the subsequences and reports
#' only those contigs that reach a finely tuned threshold of identity and coverage against de db.
#'
#'
#' @param query path to query multifasta file
#' @param database /path/to/minimap2/database/file.mmi
#' @param result name of the folder to write with the results obtained
#' @param fsize length of the subsequences to be mapped (1500 minimum)
#' @param wsize length of the sliding window to generate subsequences
#' @param rm.interm LOGICAL remove intermediate files?
#'
#' @example
#'
#' @seealso sum.info() used to merge outputs of repsearch() and minidist()
#' from the same dataset
#'
#'

minidist <- function(

  query,         # Query multifasta file
  database,      # path/to/minimap2/database/file.mmi
  result,        # Name of the file to save the results
  fsize=500,    # Size of fragments to be generated from contigs
  wsize=500,     # Sliding window size
  threads=3,
  rm.interm = T  # Remove intermdiate files (PAF tables and sequences)?

) {

  # Dependencies

  library(Biostrings, warn.conflicts = FALSE, quietly=T)
  library(seqinr, warn.conflicts = FALSE, quietly=T)
  library(dplyr)

  # Change contigs names

  #Rename contigs and move files

  fst  <- read.fasta(query)
  ctr  <- 0
  nm.s <- character(0)

  for (x in 1:length(fst)) {

    ctr  <- ctr+1
    nm   <- paste("Contig-100_",ctr, sep = "")
    nm.s <- c(nm.s, nm)

  }

  write.fasta(fst, nm.s, file.out = "tmp.fa")

  rest1  <- result
  tlsr   <- strsplit(result, split='/')[[1]]
  result <- tlsr[length(tlsr)]


  # Define internal function spliter

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

  q1 <- readDNAStringSet("tmp.fa")
  q  <- q1[width(q1)>500]

  l <- length(q)


  qnames       <- c()
  similarities <- list()
  ctg <- c()
  sbj <- c()
  system(paste0('mkdir tmp'))


  l <- length(q)
  qns <- names(q)

  dnm <- tibble()

  #Aca hacer multithread

  for (x in 1:l) {

    fout <- tempfile(tmpdir='./tmp',fileext='.fasta')
    pout <- gsub('.fasta','.paf',fout)


    s  <- as.list(as.character(q[x]))
    n  <- names(s)
    n2 <- strsplit(n,' ')[[1]][1]
    v  <- lapply(s,s2c)[[1]]

    qnames <- c(qnames,n2)

    spl <- spliter(w=fsize,s=wsize,v=v)
    len <- length(spl)
    nms <- paste(n2,1:len,sep='~')

    write.fasta(spl,names=nms,file=fout)

    vct  <- c(fout, pout,n)
    dnm[x,1:3]  <- vct

    cnc <- paste0('cat ',fout,' >> tmp/todos.fasta')
    system(cnc)
  }

  fouts <- dnm$V1
  pouts <- dnm$V2
  nouts <- dnm$V3

  cmd <- paste0('minimap2 -x asm10 ',database,' ./tmp/todos.fasta -2 -t ',threads,' > ./tmp/todos')
  system.time(system(cmd))

  # Load results: PAF table

  if (file.info('tmp/todos')$size>0) {

    tab.all <- read.csv('tmp/todos',sep='\t',header=F)

    colnames(tab.all) <- c('qid','qlen','qst','qend','strand','sid','slen','sst','send','match','len','qual')

    lnm <- length(nouts)

    for (x in 1:lnm) {

      # Parse PAF table
      tnc <- nouts[x]
      tn  <- paste0(tnc,'~')
      tab <- tab.all[grep(tn, tab.all$qid),]

      pid  <- (tab$match/tab$len)*100
      qcov <- (tab$len/tab$qlen)
      S    <- (pid*qcov)
     if (dim(tab)[1]==0) {

       similarities[[x]] <- 0

    } else {
      similarities[[x]] <- S

      dtf <- as.data.frame(table(tab$sid))
      csx <- dtf[order(dtf$Freq, decreasing = TRUE),]
      dbh <- csx$Var1[1]
      tq  <- tab$qid
      dbq <- tq[1]

      ctg[x]  <- as.character(as.vector(dbq))
      sbj[x]	<- as.character(as.vector(dbh))
    }
  }

  } else {

    similarities[[x]] <- rep(0,len)

  }


    # Remove intermediate files?

    if (rm.interm == T) {

      system(paste('rm -rf ./tmp'))
    }


  names(similarities) <- qnames
  rslt <- paste0(result,".tsv")
  Rslt <- paste0(result,".RDS")

#  saveRDS(similarities, file = Rslt )


# Filtering mapped contigs id.qcov > 65 for plasmids

 mean.vect <- rapply(similarities, mean)
 sim       <- similarities[mean.vect>=65]
 ctg2      <- as.vector(na.omit(names(sim)))

 saveRDS(sim, file = Rslt )

#Saving results in output table



 plsdb_table <- readRDS("Data/plsdb_table.RDS")
 plsdb_names <- plsdb_table$plasmid_name

    cntg2   <- rep(NA, length(ctg))
    lnth2   <- rep(NA, length(ctg))
    ht.vct2 <- rep(NA, length(ctg))



    for (i in 1:length(ctg)){

      vn1     <- strsplit(ctg,split="~")[[i]][1]
      vn      <- vn1 %in% ctg2

      if (vn == TRUE) {
      cntg2[i]   <- vn1
          qry    <- q[grep(vn1, qns, fixed=TRUE)]
      lnth2[i]   <- width(qry)[1]
      ht.vct2[i] <- grep(sbj[i], plsdb_names, value = TRUE)
      } else {
        cntg2[i]   <- NA
        lnth2[i]   <- NA
        ht.vct2[i] <- NA
      }
   }

      cntg   <- as.vector(na.omit(cntg2))
      lnth   <- as.vector(na.omit(lnth2))
      ht.vct <- as.vector(na.omit(ht.vct2))

      #Save results in table

      dfn <- data.frame()
  for (i in 1:length(ht.vct)) {

     ht <- ht.vct[i]

     dfg <- subset(plsdb_table, plsdb_table$plasmid_name == ht, select = c(plasmid_name, plasmid_length))
     dfn <- rbind(dfn,dfg)

  }

  pls_db_nm  <- as.vector(dfn$plasmid_name)
  pls_db_len <- as.vector(dfn$plasmid_length)


  dfv        <- data.frame(cntg, lnth, pls_db_nm, pls_db_len)
  names(dfv) <- c("contigs", "length", "plsdb_name", "plsdb_length")

 #Agregar sep = "\t" asi queda con un tab y puedo mejorar el parseo en sum.info
  write.table(dfv, file = rslt, quote = FALSE)
  write(cntg, file = "list.hitsdb.txt")

  #Recover contigs

  fst <- paste0("plasmid_cntgs.fasta")

  mkbldb.cntgs <- paste("makeblastdb -in tmp.fa -dbtype nucl -out ",query,".db ", "-parse_seqids -hash_index", sep ="")

  system(mkbldb.cntgs, wait=T, ignore.stdout = T)

  blcmd1 <- paste("blastdbcmd -db ",query,".db -entry_batch list.hitsdb.txt > ",fst, sep= "")

  system(blcmd1, wait = T)

  ff   <- paste("mkdir ",rest1, sep='')
  mvfr <- paste("mv", fst, rslt, Rslt, rest1, sep= " ")

  system(ff, ignore.stdout = T, ignore.stderr = T)
  system(mvfr, ignore.stdout = T, ignore.stderr = T)

  if (rm.interm == T) {

    system(paste('rm -rf ',query,".db* ","list.hitsdb.txt", "tmp.fa", sep = ""))
  }


}

