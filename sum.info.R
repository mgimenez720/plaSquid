
#' sum.info is a function that merges the results of minidist and repsearch
#' @param minidist.out the result folder of minidist() run
#' @param repsearch.out the result folder of repsearch() run
#' @param out_file name of the file to be written with the merged result of both runs
#'
#' @example
#'
#' @author Matías
#'
#' @details sum.info uses the package tidyverse to merge the .csv reports of the afforementioned
#'          functions. The output of this function will not repeat the information of the contigs
#'          that were found by both functions.
#'
#'
#'
 sum.info <- function(minidist.out,   #Result folder of minidist run
                      repsearch.out,  #Result folder of repsearch run
                      out_dir,
                      query
            ) {


   library(readr)
   library(tidyverse)
   library(Biostrings)

  #Parse minidist result



   tlsr  <- strsplit(minidist.out, split='/')[[1]]
   reslt <- tlsr[length(tlsr)]
   m.dir <- paste0(minidist.out,"/",reslt,'.tsv')
   m.fas <- paste0(minidist.out,"/plasmid_cntgs.fasta")
   m.csv <- read_delim(m.dir, delim=' ')
   names(m.csv)<- c('X1', 'contigs','length','plsdb_name')
   m <- select(m.csv, contigs, length, plsdb_name)

   sdr <- readRDS(paste0(minidist.out,'/',reslt,'.RDS'))
   S.value <- rapply(sdr, mean)
   mv <-  data.frame(S.value, contigs = as.character(names(S.value)))

   m1 <- left_join(m, mv, by = c('contigs' = 'contigs'))

  #Parse repsearch result

   r.dir <- paste(repsearch.out,"/report.csv", sep = '')
   r.csv <- read_csv(file = r.dir)
   r <- select(r.csv, contigs, Rep_type, MOB_group)

  #Summary table

   m2 <- full_join(x = r, y = m1, by = "contigs")
   system(paste0('mkdir ',out_dir), ignore.stderr = TRUE)


  #Extract plasmidic contigs

  # write(hts, file=paste0(out_dir,"/contigs.txt"))

   hts <- m2$contigs

   l <- length(hts)

   cntgs <- c()

    for (i in 1:l) {

      ht <- hts[i]

      cntgs[i] <- as.numeric(strsplit(ht, split ="_")[[1]][2])


    }

      fst  <- readDNAStringSet(query)
      h    <- fst[cntgs]
      fout <- paste0(out_dir,'/all_plasmids.fasta')
      writeXStringSet(h, fout)

      cont_names <- names(h)

      m2 <- cbind(m2, cont_names)

      write.table(m2, file=paste0(out_dir,"/Plasmids_all.csv"), sep= "\t")

      return(m2)
}

