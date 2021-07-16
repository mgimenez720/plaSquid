#!/usr/bin/env nextflow
   
process Splitter {

  label 'big_mem'
  label 'big_cpus'
    
    input:
    path "spl_contigs"
     
    output:
    path "plasmid.split"

    script:         
    """
    
    Rscript $baseDir/bin/splitter.R spl_contigs
    
    """
}


process Mapping_pr {
    
    label 'big_mem'
    label 'big_cpus'

    input:
    path "plasmid.split.final"
    path "plsdb.mmi"

    output: 
    path "plasmid.split.paf"

    script:
    """
    
    minimap2 \
    -x asm5 \
    plsdb.mmi \
    plasmid.split.final \
    -t ${task.cpus} \
    > plasmid.split.paf

    """

}

process Parse_paf {
    label 'big_mem'
    label 'big_cpus'

    tag "$contig_id"
    
    input:
    path "plasmid.split.paf"
    path contigs
    path "plasmid.split.final"

    output:
    path "Minidist_result.tsv"

    script:
    """
    
    Rscript $baseDir/bin/Parse_paf.R \
    plasmid.split.paf \
    ${contigs} \
    plasmid.split.final \
    $baseDir/data/plsdb_table.RDS

    """
}

process RetrievePlasmids {
  label 'small_mem'
  label 'small_cpus'
  
   input:
   path "temp.tsv"
   path contigs
   
   output:
   path "Minidist_result.tsv"
   

   script:
   """
   Rscript $baseDir/bin/Retrieve_plasmids.R temp.tsv ${contigs}
   
   """
}

process MinidistOut {
  
  label 'small_cpus'

  publishDir "$params.outdir/", mode: "copy"

  input:
  path "Minidist_result.tsv"
  path contigs

  output:
  path "Result.fasta"
  path "Result.tsv"
  
  script:
  """
 
  Rscript $baseDir/bin/sum.minidist.R Minidist_result.tsv ${contigs}

  """

}

process RepsearchOut {

  label 'small_cpus'

  publishDir "$params.outdir/", mode: "copy"

  input:
  path contigs
  path "Plasmid_Report.tsv"
  

  output:
  path "Result.tsv"
  path "Result.fasta"

  script:
  """
  repsearch_sum.R Plasmid_Report.tsv ${contigs}

  """

}

process Renamecntgs {
  label 'big_mem'
  label 'small_cpus'

   input:
   path contigs
   
   output:
   path "assembly.fa"

   script:
   """
   
   awk '/^>/{print ">contig-" ++i; next}{print}' ${contigs} > assembly.fa 
   
   """
}

process AnnotContigs {
label 'big_mem'
label 'big_cpus'  

   input:
   path "assembly.fa"

   output:
   path "prots.faa"

   script:
   """

   prodigal -p meta -a prots.faa -i assembly.fa
  
   """
}

process RepSearch {
label 'big_mem'
label 'big_cpus'

   input:
   path "prots.faa" 

   output:
   path "ProtsvsRep.tsv"

   script:
   """

   hmmsearch --cut_ga --cpu ${task.cpus} -o log --domtblout ProtsvsRep.tsv $baseDir/data/All_Rep.hmm prots.faa 
   

   """
}

process DomainArch {
label 'small_mem'
label 'big_cpus'

  input:
  path "ProtvsRep.tsv"

  output:
  path "multi_dom_RIP.tsv"
  path "single_dom_RIP.tsv"
  path "Domain_Architecture.RDS"

  script:
  """

  Rscript $baseDir/bin/Dom_Arch.R ProtvsRep.tsv

  """

} 

process FilterDom {
label 'small_mem'
label 'big_cpus'

  input:
  path "multi_dom_RIP.tsv"
  path "single_dom_RIP.tsv"
  path "Domain_Architecture.RDS"
  

  output:
  path "Rep_domains.tsv"

  script:
  """

  Rscript $baseDir/bin/Filter_RIP.R multi_dom_RIP.tsv single_dom_RIP.tsv Domain_Architecture.RDS $baseDir/data/Arq_RIP_new.RDS

  """

} 

process IncSearch {

label 'big_mem'
label 'big_cpus'

  input:
  path "prots.faa"

  output:
  path "Inc_candidates.tsv"

  script:
  """
  hmmsearch -o log --cpu ${task.cpus} --domtblout Inc_candidates.tsv $baseDir/data/All_Inc.hmm prots.faa 

  """

}

process RnaSearch {
label 'big_mem'
label 'big_cpus'

  input: 
  path "assembly.fa"

  output:
  path "RNA_candidates.tsv"

  script:
  """

  cmsearch --cpu ${task.cpus} --tblout RNA_candidates.tsv $baseDir/data/All_RNA_Inc.hmm assembly.fa 
  
  """
}

process IncClassif {

label 'small_mem'
label 'small_cpus'

  input: 
  path "Inc_candidates.tsv"
  path "RNA_candidates.tsv"

  output:
  path "Classification_table.tsv"

  script:
  """

  Rscript $baseDir/bin/Inc_classification.R Inc_candidates.tsv RNA_candidates.tsv 
  
  """
}

process IncFilter {

label 'small_mem'
label 'small_cpus'
  
  input: 
  path "Classification_table.tsv"

  output:
  path "Filtered_Classif.tsv"

  script:
  """

  Rscript $baseDir/bin/Filter_classification.R Classification_table.tsv 
  
  """
}

process MobSearch {
label 'big_mem'
label 'big_cpus'

  input:
  path "prots.faa"

  output:
  path "Mob_candidates.tsv"

  script:
  """
  hmmsearch -o log --cpu ${task.cpus} --domtblout Mob_candidates.tsv $baseDir/data/All_MOB.hmm prots.faa 

  """

}

process MobFilter {
label 'small_mem'
label 'small_cpus'

  input: 
  path "Mob_candidates.tsv"

  output:
  path "Mob_table.tsv"

  script:
  """

  Rscript $baseDir/bin/Filter_Mob.R Mob_candidates.tsv 
  
  """
}

process GeneRetrieve {

label 'big_mem'
label 'big_cpus'

  input:
  path "Filtered_classif.tsv"
  path "Mob_table.tsv"
  path "Rep_domains.tsv"
  path "assembly.fa"

  output:
  path "Plasmids_contigs.fasta"
  path "Plasmid_Report.tsv"

  script:
  """
  
  Rscript $baseDir/bin/Retrieve_RIP_plasmids.R Filtered_classif.tsv Rep_domains.tsv Mob_table.tsv assembly.fa

  """

}
process SumOutput {

label 'small_mem'
label 'big_cpus'

publishDir "$params.outdir/", mode: "copy"

  input:
  path "Minidist_result.tsv"
  path "Plasmid_Report.tsv"
  path contigs

  output:
  path "Result.fasta"
  path "Result.tsv"

  script:
  """
  
  Rscript $baseDir/bin/sum.info.R Minidist_result.tsv Plasmid_Report.tsv ${contigs} 
  
  """

}

process DownPLSDB {

label 'big_mem'
label 'big_cpus'

  output:
  path "plsdb.fna"
  
  script:
  """
  
  wget -O plsdb https://ndownloader.figshare.com/files/23582252
  unzip plsdb
  rm PLSDB_2020_06_29/plsdb.msh PLSDB_2020_06_29/plsdb.tsv PLSDB_2020_06_29/plsdb.abr
  blastdbcmd -dbtype nucl -db PLSDB_2020_06_29/plsdb.fna -entry all -out plsdb.fna
    
  """

}

process FormtPLSDB {
label 'big_mem'
label 'big_cpus'

  input:
  path "plsdb.fna"

  output:
  path "plsdb.mmi"

  script:
  """
  
  minimap2 -d plsdb.mmi plsdb.fna 
  cp plsdb.mmi $baseDir/plsdb.mmi

  """
}
