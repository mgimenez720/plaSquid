#!/usr/bin/env Rscript


#' \strong{rePsearch}: Replicon and MOB search and classification
#'
#' @description  this function looks for plasmid derived sequences in metagenomic contigs
#' through a conserved genes approach. It looks for \strong{replication initiator proteins}
#' and \strong{relaxases} domains. RIPs are also filtered through domain architecture.
#' Plasmidic contigs are extracted in a fasta file and classifyed through HMMs in replicon
#' and MOB types. Results are a table indicating plasmidic contigs and its clasification.
#'
#' @param file Fasta file with metagenomic assembled contigs.
#' @param out Name of the folder to create and write results.
#' @param rm.interm Logical if T removes intermediate files.
#'
#' @seealso sum.info() used to merge outputs of repsearch() and minidist()
#' from the same dataset
#'
#' @example
#'
#'
#' rep_models contains hmms corresponding to replication initiator proteins domains in pfam.
#' Inc_models contains hmms for classification in incompatibility groups
#' Rel_models contains hmms corresponding to relaxases domains in pfam.
#' Mob_models contains hmms for classification in MOB groups.
#' Arch.rep.list.RDS contains the list of domain architectures of replication initiator proteins.

#contigs annotation


repsearch <- function(

  file,
  out,
  rm.interm

)          {

  #dependencies

  library('seqinr', warn.conflicts = FALSE)
  library('dplyr', warn.conflicts = FALSE)
  library('tidyverse', warn.conflicts = FALSE)



  #Rename contigs and move files

  fst <- read.fasta(file)
  ctr <- 0
  nm.s <- character(0)

  for (x in 1:length(fst)) {

    ctr  <- ctr+1
    nm   <- paste("Contig-100_",ctr, sep = "")
    nm.s <- c(nm.s, nm)


  }


  write.fasta(fst, nm.s, file.out = "tmp.fa")

  s     <- strsplit(file, "/", fixed = TRUE)[[1]]
  file  <- s[length(s)]

  system(paste("mkdir ",out))

  mv1 <- paste("cp"," -r ","./Data/Rep_models ","./",out,"/Rep_models", sep='')
  mv2 <- paste("cp"," -r ","./Data/Inc_models ","./",out,"/Inc_models", sep='')
  mv3 <- paste("cp"," -r ","./Data/Rel_models ","./",out,"/Rel_models", sep='')
  mv4 <- paste("cp"," -r ","./Data/Mob_models ","./",out,"/Mob_models", sep='')
  mv5 <- paste("cp ./tmp.fa ./",out,"/",file, sep='')

  system(mv1,ignore.stderr=T,ignore.stdout=T)
  system(mv2,ignore.stderr=T,ignore.stdout=T)
  system(mv3,ignore.stderr=T,ignore.stdout=T)
  system(mv4,ignore.stderr=T,ignore.stdout=T)
  system(mv5,ignore.stderr=T,ignore.stdout=T)
  system("rm tmp.fa", ignore.stderr = T, ignore.stdout = T)
  setwd(out)

  message("Annotating contigs...")

  #Get orfs in the metagenome

  orf <- paste("prodigal -q"," -i ",file ," -a ","./",file,".faa",sep='')

  system(orf, wait = T, ignore.stdout=T)

  prot <- paste(file,".faa", sep="")
  orfs <- read.fasta(file = prot)
  annotation <- getAnnot(orfs)

  #Make blast databases

  mkbldb.faa   <- paste("makeblastdb -in ", prot," -dbtype prot -out ",prot,".db ", "-parse_seqids -hash_index", sep ="")
  mkbldb.cntgs <- paste("makeblastdb -in ", file," -dbtype nucl -out ",file,".db ", "-parse_seqids -hash_index", sep ="")

  system(mkbldb.faa, wait=T, ignore.stdout = T)
  system(mkbldb.cntgs, wait=T, ignore.stdout = T)


  #Search Replication Initiator protein domains

  message("looking for replication initiator proteins...")

  #lst.rd <- paste("cp ./../Arch.rep.list.RDS ./Arch.rep.list.RDS", sep = "")
  #system(lst.rd, wait=T, ignore.stdout = T)
  lst <- readRDS("./../Data/Arch.rep.list.RDS")

  setwd("./Rep_models")

  reps <- c("IncFII_repA", "PriCT_1", "Rep_1", "Rep_2", "Rep_3", "RepA_C", "RepA_N", "RepC", "RepL", "Rep_N", "Replicase", "Rep_trans",
            "RHH_1", "Rop", "RPA", "RP_C", "TrfA")


  Srchr <- function(out_dir,      #Path to output director
                    Mdls,         #R vector with HMM models names
                    out_lst,      #name of the file with orf headers
                    out_orfs      #name of the output fasta file with orfs
  )
  {

    system(paste("mkdir ",out_dir, sep=""))
    lst.orfs <- character(0)

    for (x in 1:length(Mdls)) {

      #Searching for orfs with hmm domains

      Mdl <- Mdls[x]
      Domain.srch <- paste('hmmsearch --cut_ga -o ./log --tblout ',out_dir,'/tablevs',Mdl,' ', Mdl,'.hmm ./../',file,'.faa' ,sep='')

      system(Domain.srch,ignore.stderr=T,ignore.stdout=T)

      #Retrieve Relaxases' orfs.

      table.hmm <- read.table(paste(out_dir,"/tablevs",Mdl, sep =""), header = FALSE, sep = "", blank.lines.skip = TRUE, skipNul = TRUE,
                              col.names= c("contigs","taccession","queryname", "qaccession","Evalue","score","bias","domEvalue","domscore",
                                           "dombias","exp","reg","clu","ov","env","dom","rep","inc"))

      if ((is.data.frame(table.hmm) && nrow(table.hmm)==0)==TRUE){

        remove(table.hmm)

        next

      }

      name1    <- as.vector((table.hmm$contigs), mode ="any")
      lst.orfs <- c(lst.orfs, name1)

    }
    #Concatenate hits lists and retrieve orfs

    setwd(out_dir)

    write(lst.orfs, file = paste(out_lst, sep = ""),sep = "\n")

    blcmd <- paste("blastdbcmd -db ./../",prot,".db -entry_batch ",out_lst," > ",out_orfs, sep= "")

    system(blcmd, wait=T)



  }

  Srchr(out_dir = "./../rep_search",
        Mdls = reps,
        out_lst = "list.hits.Rep",
        out_orfs = "hits.Rep.faa")

  if (file.size("hits.Rep.faa") == 0){

    warning("No replication initiator proteins found")

  } else {

    #Search other domains in RIP orfs.

    hmm.srch <- paste("hmmsearch --cut_ga --noali -o.log --domtblout ./tableallvsRep ./../Rep_models/all.hmm ./hits.Rep.faa", sep="")
    system(hmm.srch, wait=T, ignore.stdout = T)

    all.hmm.tab <- read.table("tableallvsRep", header = FALSE, fill=TRUE,
                              col.names= c("contigs","taccession","tlen","queryname",
                                           "qaccession","qlen","Evalue","score","bias",
                                           "num","of","CEvalue","iEvalue","domscore","dombias",
                                           "hmmfrom","hmmto","alifrom","alito","envfrom",
                                           "envto","acc"), sep = "")

    #Parse multiple domains or single domain proteins.
    #Usar herramientas de dplyr para encarar esta parte!!!!! Hay errores!

    all.hmm.tab %>% group_by("contigs")

    multi.domain.tab1  <- all.hmm.tab[duplicated(all.hmm.tab$contigs) | duplicated(all.hmm.tab$contigs, fromLast=TRUE), ]
    multi.domain.tab   <- multi.domain.tab1[unique(multi.domain.tab1$queryname),]
    mdt                <- multi.domain.tab[with(multi.domain.tab, order(contigs, envfrom)), ]
    contigs            <- unique(mdt$contigs)


    unique.domain.tab1 <- anti_join(all.hmm.tab, multi.domain.tab, by='contigs')
    unique.domain.tab  <- distinct(unique.domain.tab1, unique.domain.tab1$contigs, .keep_all = TRUE)


    #Computing domain architecture of multi domain proteins

    m <- character(0)

    for (x in 1:length(contigs)) {

      cntg    <- contigs[x]

      contig  <- subset.data.frame(mdt, mdt$contigs==cntg, select = c(contigs, queryname, score, envfrom, envto))
      contig1 <- contig[order(contig$envfrom),]
      contig2 <- subset(contig1, !duplicated(subset(contig1, select=c(contigs, queryname, envfrom))))

      e1 <- contig2[1,5]
      s2 <- contig2[2,4]

      if (is.na(s2)==TRUE){
        dom1a <- as.character(contig2[1,2])
        dom1  <- as.vector(c(dom1a,"no_2nd_domain"))
      } else {
        if (e1<s2){
          dom1a <- as.character(contig2[1,2])
          dom1  <- as.vector(c(dom1a, "not_over"))
        } else {
          bs1 <- contig2[1,3]
          bs2 <- contig2[2,3]
          if (bs1<bs2){
            dom1a <- as.character(contig2[2,2])
            dom1  <- as.vector(c(dom1a,"overlapped"))
          } else {
            dom1a <- as.character(contig2[1,2])
            dom1  <- as.vector(c(dom1a,"overlapped"))
          }
        }
      }

      #Second overlap

      e2 <- contig2[2,5]
      s3 <- contig2[3,4]

      if (is.na(s3)==TRUE){
        dom2a <- as.character(contig2[2,2])
        dom2  <- as.vector(c(dom2a,"no_3rd_domain"))
      } else {
        if (e2<s3){
          dom2a <- as.character(contig2[2,2])
          dom2  <- as.vector(c(dom2a, "not_over"))
        } else {
          bs2 <- contig2[2,3]
          bs3 <- contig2[3,3]
          if (bs2<bs3){
            dom2a <- as.character(contig2[3,2])
            dom2  <- as.vector(c(dom2a,"overlapped"))
          } else {
            dom2a <- as.character(contig2[2,2])
            dom2  <- as.vector(c(dom2a,"overlapped"))
          }
        }
      }

      #Third overlap

      e3 <- contig2[3,5]
      s4 <- contig2[4,4]

      if (is.na(s4)==TRUE){
        dom3a <- as.character(contig2[3,2])
        dom3  <- as.vector(c(dom3a,"no_4th_domain"))
      } else {
        if (e3<s4){
          dom3a <- as.character(contig2[3,2])
          dom3  <- as.vector(c(dom3a,"not_over"))
        } else {
          bs3 <- contig2[3,3]
          bs4 <- contig2[4,3]
          if (bs3<bs4){
            dom3a <- as.character(contig2[4,2])
            dom3  <- as.vector(c(dom3a,"overlapped"))
          } else {
            dom3a <- as.character(contig2[3,2])
            dom3  <- c(dom3a,"overlapped")
          }
        }
      }

      #Fourth and last overlap

      e4<-contig2[4,5]
      s5<-contig2[5,4]

      if (is.na(s5)==TRUE){
        dom4a <- as.character(contig2[4,2])
        dom4  <- as.vector(c(dom4a,"no_5th_domain"))
      } else {
        if (e4<s5){
          dom4a <- as.character(contig2[4,2])
          dom4  <- as.vector(c(dom4a,"not_over"))
        } else {
          bs4 <- contig2[4,3]
          bs5 <- contig2[5,3]
          if (bs4<bs5){
            dom4a <- as.character(contig2[5,2])
            dom4  <- as.vector(c(dom4a,"overlapped"))
          } else {
            dom4a <- as.character(contig2[4,2])
            dom4  <- as.vector(c(dom4a,"overlapped"))
          }
        }
      }


      Arq <- as.vector(c(dom1,dom2,dom3,dom4))

      #Removing NAs from domain architecture vector

      Arq <- Arq[!is.na(Arq)]


      sapply(lst,identical,Arq)->comparison


      if ("TRUE" %in% comparison==TRUE){
        h <- as.vector(cntg)
        m <- append(m, h)


      }

    }


    write(m, file = "result_double_all")
    dd <- unique(scan(file = "result_double_all", what = "character", sep = "\n"))
    write(dd, file="list.multi.doms")

    write.table(unique.domain.tab, file="Single_tab_all", sep = " ", row.names = TRUE, col.names = TRUE)



    #Filtering single-domain RIPs by bit-score and length.


    lsds  <- c("IncFII_repA","RepA_C", "RepA_N", "RepC", "Replicase", "Rop", "RPA", "RP-C", "TrfA", "IncP9", "IncAC")
    lsdvs <- c("10.0","45", "38", "76", "76", "77.1", "67.8", "45", "87", "236.4", "655.7")

    names(lsdvs) <- lsds

    sd1 <- character(0)

    for (x in 1:length(lsds)) {

      lsd   <- lsds[x]

      s_dom <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == lsd)
      write.table(s_dom, file = paste0("single.tab_",lsd))

      if (is.na(file.size(paste0("single.tab_",lsd)))) {

        next
      } else {

        single.domain   <- read.table(paste0("single.tab_",lsd), header = TRUE, fill=TRUE, sep = "")
        tab.single.hits <- single.domain[(single.domain$score > as.numeric(lsdvs[lsd])) ,]
        sd.list         <- as.vector(tab.single.hits$contigs)

        write(sd.list, file = paste0("list.hits.",lsd), ncolumns = 1, sep = "" )

        sd1 <- c(sd1,sd.list)
      }

    }


    Sdm <- c("PriCT_1","Rep_1","Rep_2","Rep_3","RepL","Rep_trans")

    sdom        <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "PriCT_1")
    PriCT_1     <- sdom[(sdom$score > 49 & sdom$tlen < 500 & sdom$tlen > 420) ,]
    PriCT_1_sdl <- as.vector(PriCT_1$contigs)

    sdom      <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "Rep_1")
    Rep_1A    <- sdom[(sdom$score > 37 & sdom$tlen > 130) ,]
    Rep_1B    <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_1     <- rbind(Rep_1A, Rep_1B)
    Rep_1_sdl <- as.vector(Rep_1$contigs)

    sdom      <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "Rep_2")
    Rep_2_sdl <- as.vector(sdom$contigs)

    sdom      <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "Rep_3")
    Rep_3_sdl <- as.vector(sdom$contigs)

    sdom     <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "RepL")
    RepL     <- sdom[(sdom$score > 85 & sdom$tlen > 90) ,]
    RepL_sdl <- as.vector(PriCT_1$contigs)

    sdom          <- subset.data.frame(unique.domain.tab, unique.domain.tab$queryname == "Rep_trans")
    Rep_trans     <- sdom[(sdom$score > 27 & sdom$tlen < 130) ,]
    Rep_trans_sdl <- as.vector(Rep_trans$contigs)



    sd2 <- c(PriCT_1_sdl, Rep_1_sdl, Rep_2_sdl, Rep_3_sdl, RepL_sdl, Rep_trans_sdl)
    sd  <- c(sd1, sd2)
    hit.contigs <- c(dd, sd)

    write(hit.contigs, file = "list.rep.hits")

    #Retrieve contigs encoding RIP

    Retriever <- function( m,      #path to file of ORFs list saved as R vector
                           fldb,   #path to blastn database with contigs to retrieve
                           cntnm,  #file that will contain the contigs headers to be retrieved
                           fst     #path to Fasta file to be saved
    )
    {

      cntgs.tab <- read.table(m, header = FALSE, col.names = "headers")
      df        <- separate(cntgs.tab, headers, c('head','cont','orf'),sep = '_',remove = FALSE, extra = "drop")
      df2       <- unite(df, contigs, head, cont, sep = '_', remove = FALSE)

      write(unique(as.vector(df2$contigs)), file = cntnm)

      blcmd1 <- paste("blastdbcmd -db ",fldb," -entry_batch ",cntnm," > ",fst, sep= "")

      system(blcmd1, wait=T)


    }

    Retriever(m = "list.rep.hits", fldb = paste0("./../",file,".db"), cntnm = "headers.rep", fst = "contigs.rep.fasta")




    #Retrieving filtered RIP

    bcm <- paste("blastdbcmd -db ./../",prot,".db"," -dbtype prot -entry_batch list.rep.hits -out rep.hits.faa", sep="")

    system(bcm,ignore.stderr=T,ignore.stdout=T, wait=T)

    #RIP classification with Inc_models

    system(paste("mkdir ./../rep_classification", sep=""))

   # Incs <- c("IncA-C", "IncB-O_IncK_IncI","IncFIA","IncFIB","IncFII_IncFIC_IncZ","IncG-U","IncHI2","IncHIA","IncHIB",
    #          "IncL-M", "IncN","IncP-1", "IncP-2", "IncP-7", "IncP-9", "IncQ","IncR","IncT","IncW","IncX","Ref_rep10",
     #         "Ref_rep11", "Ref_rep12", "Ref_rep13", "Ref_rep14", "Ref_rep15", "Ref_rep16", "Ref_rep17", "Ref_rep18",
     #         "Ref_rep19","Ref_rep20","Ref_rep21","Ref_rep22","Ref_rep23","Ref_rep24","Ref_rep1","Ref_rep2","Ref_rep3",
    #        "Ref_rep4","Ref_rep5","Ref_rep6","Ref_rep7","Ref_rep8","Ref_rep9")


Incs <- c("Inc10", "Inc13", "Inc18", "Inc7", "Inc9", "Inc1", "IncB-O-K-I", "IncFII_RepE", "IncFI_RepE", "IncLM", "IncP2",
          "IncP9", "IncQ", "IncT", "IncX", "Inc11", "Inc14", "Inc4", "Inc8", "IncAC", "IncFII_RepA1", "IncFI_RepB",
          "IncH", "IncGU", "IncN", "IncP7", "IncP", "IncR", "IncW", "IncZ")

    #scores.Inc <- c(710.8, 569.5, 606.1, 731.0, 588.8, 849.8, 855.2, 650.8, 683.6, 592.0, 317.0, 614.9, 674.4, 679.3,
     #               256.1, 577.6, 658.1, 675.5, 563.5, 390.0, 328.9, 598.7, 1233.1, 415.4, 589.2, 742.2, 533.3, 806.1,
      #              389.6, 548.9, 361.2, 530.0, 561.7, 789.7, 681.2, 312.7, 1196.6, 1167.0, 611.4, 414.3, 481.3, 567.2,
       #             672.7, 619.2)

scores.Inc <- c(307.7, 248.6, 453.4, 138.3, 302.6, 301.6, 573.9, 825.5, 517.8, 466.2, 392.1, 236.4,
                205.7, 575.2, 416.6, 164.3, 406.4, 384.6, 657.4, 655.7, 664.9, 438.0, 173.8, 654.2, 257.0,
                674.8, 571.7, 247.3, 269.4, 515.0)


    names(scores.Inc) <- Incs
                  hre <- list()

    for (x in 1:length(Incs)) {

      Inc <- Incs[x]

      hmm.srch.1 <- paste("hmmsearch --noali -o log --domtblout ./../rep_classification/hit.repvs",Inc,
                        " ./../Inc_models/",Inc,".hmm ./rep.hits.faa", sep="")

      system(hmm.srch.1, wait=T, ignore.stdout = T)

      tabla <- read.table(paste("./../rep_classification/hit.repvs",Inc, sep = ""), header = FALSE, fill=TRUE,
                          col.names= c("contigs","taccession","tlen","queryname", "qaccession","qlen","Evalue","score",
                                       "bias","num","of","CEvalue","iEvalue","domscore","dombias","hmmfrom","hmmto",
                                       "alifrom","alito","envfrom","envto","acc", "DoT"))


      hits.rep <- subset.data.frame(tabla, tabla$score>=unname(scores.Inc[Inc]), select = c(contigs, queryname))


      if ((is.data.frame(hits.rep) && nrow(hits.rep)==0)==TRUE){

        next

      }  else {

        hits.rep.exp   <- separate(hits.rep, contigs, c('head','cont','orf','orf2'),sep = '_')
        hits.rep.exp2  <- unite(hits.rep.exp, contigs, head, cont, sep = '_', remove = FALSE)
        hits.rep.final <- unite(hits.rep.exp2, newcontigs, head, cont, queryname, sep = '_')

        #####CAMBIAR EL NOMBRE de la tabla y hacerlo salida!

        write.table(hits.rep.final, file = paste("./../rep_classification/table.hits.",Inc, sep = ""),
                    append = FALSE, col.names = TRUE)

        hre[[x]] <- tibble(contigs = hits.rep.exp2$contigs, Rep_Type = as.vector(hits.rep.exp$queryname))

        #write.table(hre, file = paste0("./../table.hits.",Inc),append = FALSE, col.names = c("contigs", "Inc_group"))

      }
    }
        if (identical(hre,list())) {

          warning('new RIPs found')


          } else {

        cdf <- bind_rows(hre)

        colnames(cdf)  <- c("contigs", "Rep_type")

        write.table(cdf, file = paste0("./../Inc_table"), append = FALSE, col.names = c("contigs", "Rep_type"))
        }

    # Change headers with Inc groups classification



    Classifyer <- function(fst,     #Path to fasta file to change headers
                           grps,    #Vector with the names of the groups to add to the headers
                           nhr,     #Path to table with contigs classifyed in each group
                           fstout   #Path to fasta file with new headers
    )
    {

      f  <- read.fasta(file = fst, forceDNAtolower = FALSE)
      n  <- names(f)
      n2 <- paste0(n,' ')


      for (x in 1:length(grps)) {

        grp <- grps[x]

        if (is.na(file.size(paste0(nhr,grp)))) {

          next
        }
        chngr    <- read.table(paste0(nhr,grp), header=TRUE)
        Grpcntgs <- as.vector(chngr$contigs)

        for (j in 1:length(Grpcntgs)) {


          Grpcnt <- Grpcntgs[j]
          n2     <- gsub(paste0(Grpcnt," "),paste0(Grpcnt,"__",grp), n2 , fixed= TRUE)

        }
      }

      write.fasta(f, n2, fstout)
    }

    Classifyer(fst    = "contigs.rep.fasta",
               grps   = Incs,
               nhr    = "./../rep_classification/table.hits.",
               fstout = "./../Rep.contigs.fasta")

  }

  #Relaxases search

  message("looking for relaxases...")

  setwd("./../Rel_models")
  rel.hmms <- c("MobA_MobL","MobC","Mobilization_B","Mob_Pre","TraA","TraC","TraI", "TrwC")


  Srchr(out_dir = "./../Rel_search",
        Mdls = rel.hmms,
        out_lst = "list.hits.Rel",
        out_orfs = "hits.Rel.faa")


  if (file.size("hits.Rel.faa") == 0){

    warning("No relaxases found, finishing pipeline")

  } else {



    #Classify Relaxases through hmm approach

    system(paste("mkdir ./../rel_classification", sep=""))


    Mobs       <- c("MOBC","MOBF","MOBP","MOBHEN","MOBH","MOBQ","MOBV")
    scores.MOB <- c(96.6, 95.7, 46.7, 276.0, 203.7, 63.3, 83.5)
    names(scores.MOB) <- Mobs

    for (x in 1:length(Mobs)) {

      MOB <- Mobs[x]
      hmm.srch.2 <- paste("hmmsearch --noali -o log --domtblout ./../rel_classification/hit.relvs",MOB,
                          " ./../Mob_models/",MOB,".hmm ./hits.Rel.faa", sep="")

      system(hmm.srch.2, wait=T, ignore.stdout = T)

      tabla <- read.table(paste("./../rel_classification/hit.relvs",MOB, sep = ""), header = FALSE, fill=TRUE,
                          col.names= c("contigs","taccession","tlen","queryname", "qaccession","qlen","Evalue",
                                       "score","bias","num","of","CEvalue","iEvalue","domscore","dombias","hmmfrom",
                                       "hmmto","alifrom","alito","envfrom","envto","acc", "DoT"))

      hits.Rel <- subset.data.frame(tabla, tabla$score>=unname(scores.MOB[MOB]), select = c(contigs, queryname))


      if ((is.data.frame(hits.Rel) && nrow(hits.Rel)==0)==TRUE){

        next

      }  else {

        hits.Rel.exp   <- separate(hits.Rel, contigs, c('head','cont','orf'),sep = '_', extra = "drop")
        hits.Rel.exp2  <- unite(hits.Rel.exp, contigs, head, cont, sep = '_', remove = FALSE)
        hits.Rel.final <- unite(hits.Rel.exp2, newcontigs, head, queryname, sep = '_')

        write.table(hits.Rel.final, file = paste("./../rel_classification/table.hits.",MOB, sep = ""),
                    append = FALSE, col.names = TRUE)


      }
    }

    # Retrieve contigs containing relaxases

    Retriever(   m = "list.hits.Rel",
                 fldb = paste0("./../",file,".db"),
                 cntnm = "headers.Rel",
                 fst = "contigs.rel.fasta")




    # Change headers with MOB groups classification


    Classifyer(fst = "contigs.rel.fasta",
               grps = Mobs,
               nhr = "./../rel_classification/table.hits.",
               fstout = "./../Rel.contigs.fasta")

  }

  # Remove intermediate files

  setwd("./../")

  if (rm.interm == T) {

    system(paste0('rm -rf *_search *_models *_classification ', file,'* *.RDS'))
  }

  #Generate reports


  n.rel <- as.vector(names(read.fasta(file = 'Rel.contigs.fasta')))
  n.rep <- as.vector(names(read.fasta(file = 'Rep.contigs.fasta')))

  contigs   <- character(0)
  MOB_group <- character(0)

  for (i in 1:length(n.rel)){

    nrl <- n.rel[i]

    contigs1      <- strsplit(nrl, split = "__", fixed = TRUE)[[1]][1]
    contigs[i]    <-  strsplit(contigs1, split = "|", fixed = TRUE)[[1]][2]
    MOB_group[i]  <- strsplit(nrl, split = "__", fixed = TRUE)[[1]][2]

  }

  Rel_df <- as.data.frame(cbind(contigs, MOB_group))
  write.table(Rel_df, 'report.Rel.txt')

  contigs  <- character(0)
  Rep_type <- character(0)

  for (i in 1:length(n.rep)){

    nrp <- n.rep[i]

    contigs1    <- strsplit(nrp, split = "__", fixed = TRUE)[[1]][1]
    contigs[i]  <- strsplit(contigs1, split = "|", fixed = TRUE)[[1]][2]
    Rep_type[i] <- strsplit(nrp, split = "__", fixed = TRUE)[[1]][2]

  }

  Rep_df <- as.data.frame(cbind(contigs, Rep_type))

  write.table(Rep_df, 'report.Rep.txt')

  m <- full_join(x = Rep_df, y = Rel_df, by = "contigs")

  write.csv(m, file='report.csv')

  setwd("./../")

  return(m)




}


