#!/usr/bin/env nextflow
   
process Splitter {
    
    input:
    path contigs
     
    output:
    file "Contig*"

    script:         
    """
    Rscript $baseDir/bin/script1.R ${contigs}
    """
}


process Mapping_pr {
    tag "$contig_id"

    input:
    tuple val(contig_id), file("${contig_id}.fst")

    output: 
    tuple val(contig_id), file("${contig_id}.fst"), file("${contig_id}.paf")

    script:
    """
    
    minimap2 \
    -x asm5 \
    /mnt/4tb/home/mgimenez/Matias/Metagenomas/databases/plsdb_2020.mmi \
    ${contig_id}.fst \
    > ${contig_id}.paf

    
    """

}

process Parse_paf {
    tag "$contig_id"
    
    input:
    tuple val(contig_id), file("${contig_id}.fst"), file("${contig_id}.paf")
    

    output:
    file "*.tsv"

    script:
    """
    
    Rscript $baseDir/bin/Parse_paf.R ${contig_id}.paf ${contig_id}.fst $baseDir/data/plsdb_table.RDS

    """
}


process ConcTab {

   input:
   path "*.tsv"   

   output:
   path "temp.tsv"

   script:
   """
   cat *.tsv > temp.tsv

   """
    
   
}

process RetrievePlasmids {

   input:
   path "temp.tsv"
   path contigs
   
   output:
   path "Minidist_result.tsv"
   path "Minidist_contigs.fasta"

   script:
   """
   Rscript $baseDir/bin/Retrieve_plasmids.R temp.tsv ${contigs}
   
   """
}

process AnnotContigs {

   input:
   path contigs

   output:
   path "prots.faa"
   path "assembly.fa" 

   script:
   """
   
   awk '/^>/{print ">contig-" ++i; next}{print}' ${contigs} > assembly.fa 

   prodigal -a prots.faa -i assembly.fa
  
   """
}

process RepSearch {

   input:
   path "prots.faa"
   path "assembly.fa" 

   output:
   path "ProtsvsRep.tsv"

   script:
   """

   hmmsearch --cut_ga -o log --domtblout ProtsvsRep.tsv $baseDir/data/All_Rep.hmm prots.faa 
   

   """
}

process DomainArch {

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

  input:
  path "prots.faa"
  path "assembly.fa"

  output:
  path "Inc_candidates.tsv"

  script:
  """
  hmmsearch -o log --domtblout Inc_candidates.tsv $baseDir/data/All_Inc.hmm prots.faa 

  """

}

process RnaSearch {

  input: 
  path "prots.faa"
  path "assembly.fa"

  output:
  path "RNA_candidates.tsv"

  script:
  """

  cmsearch --tblout RNA_candidates.tsv $baseDir/data/All_RNA_Inc.hmm assembly.fa 
  
  """
}

process IncClassif {

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

  input:
  path "prots.faa"
  path "assembly.fa"

  output:
  path "Mob_candidates.tsv"

  script:
  """
  hmmsearch -o log --domtblout Mob_candidates.tsv $baseDir/data/All_MOB.hmm prots.faa 

  """

}

process MobFilter {

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
  
  input:
  path "Filtered_classif.tsv"
  path "Mob_table.tsv"
  path "Rep_domains.tsv"
  path "prots.faa"
  path "assembly.fa"

  output:
  path "Plasmids_contigs.fasta"
  path "Plasmid_Report.tsv"

  script:
  """
  
  Rscript $baseDir/bin/Retrieve_RIP_plasmids.R Filtered_classif.tsv Rep_domains.tsv Mob_table.tsv assembly.fa

  """

}

process DownPLSDB {
 
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
  
  input:
  path "plsdb.fna"

  output:
  path "plsdb.mmi"

  script:
  """
  
  minimap2 -d plsdb.mmi plsdb.fna 
  
  """

  



}