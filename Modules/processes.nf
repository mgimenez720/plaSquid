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
    
    splitter.R spl_contigs
    
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
    --secondary=no \
    -2 \
    plsdb.mmi \
    plasmid.split.final \
    -t ${task.cpus} \
    > plasmid.split.paf

    """

}

process Parse_paf {
    label 'big_mem'
    label 'big_cpus'
    
    input:
    path "plasmid.split.paf"
    path contigs
    path "plasmid.split.final"

    output:
    path "Minidist_result.tsv"

    script:
    """
    
    Parse_paf.R \
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
   Retrieve_plasmids.R temp.tsv ${contigs}
   
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
 
  sum.minidist.R Minidist_result.tsv ${contigs}

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

  Dom_Arch.R ProtvsRep.tsv

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

  Filter_RIP.R multi_dom_RIP.tsv single_dom_RIP.tsv Domain_Architecture.RDS $baseDir/data/Arq_RIP_new.RDS

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

  Inc_classification.R Inc_candidates.tsv RNA_candidates.tsv 
  
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

  Filter_classification.R Classification_table.tsv 
  
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

  Filter_Mob.R Mob_candidates.tsv 
  
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
  
  Retrieve_RIP_plasmids.R Filtered_classif.tsv Rep_domains.tsv Mob_table.tsv assembly.fa

  """

}

process RipExtract {

label 'big_mem' 
label 'big_cpus'

publishDir "$params.outdir/", mode: "copy"

  input:
  path "prots.faa"
  path "Filtered_Classif.tsv"
  path "Rep_domains.tsv"

  output:
  path "RIP_seqs.faa"

  script:
  """
  
  RIP_extraction.R prots.faa Filtered_Classif.tsv Rep_domains.tsv
  
  """

}

process MobExtract {

label 'small_mem' 
label 'small_cpus'

publishDir "$params.outdir/", mode: "copy"

  input:
  path "prots.faa"
  path "MOB_table.tsv"

  output:
  path "MOB_seqs.faa"

  script:
  """
  
  MOB_extraction.R prots.faa MOB_table.tsv
  
  """

}

process SumOutput {

label 'small_mem'
label 'big_cpus'

  input:
  path "Minidist_result.tsv"
  path "Plasmid_Report.tsv"
  path contigs

  output:
  path "Result.fasta"
  path "Result.tsv"

  script:
  """
  
  sum.info.R Minidist_result.tsv Plasmid_Report.tsv ${contigs} 
  
  """

}

process FPCorrection {

label 'big_mem'
label 'small_cpus'

input:
path "Result.fasta"
path "Result.tsv"

output:
path "Chr_eval.fasta"

script:

"""

Chr_eval.R Result.fasta Result.tsv


"""

}

process FinalOut { 

  publishDir "$params.outdir/", mode: "copy"
  
  input:  
  path "Minidist_result.tsv"
  path "Plasmid_Report.tsv"
  path contigs
  
  output:
  path "plaSquid_result.fasta"
  path "plaSquid_result.tsv"
  

  script:
  '''
   
  sum.final.R Minidist_result.tsv Plasmid_Report.tsv ${contigs}

  ''' 
}

process ChrDetection {

label 'big_mem'
label 'small_cpus'

input:
path "Chr_eval.fasta"

output:
path "Bact_120.tsv"

shell:

'''
var=$(wc -c Chr_eval.fasta | cut -f1 -d " ")

if  (($var > 20000))
then
      
      prodigal -a Chr_eval.faa -i Chr_eval.fasta

      hmmsearch -o log --cut_ga --cpu !{task.cpus} --domtblout Bact_120.tsv !{baseDir}/data/Bact_120.hmm Chr_eval.faa

else

    echo "query_name" "taccession"  "tlen" "Marker_gene" "qaccession"  "qlen" "Evalue"  "score"  "bias" "num" "of"  "CEvalue"  "iEvalue"  \
    "domscore" "dombias" "hmmfrom" "hmmto" "alifrom" "alito" "envfrom" "envto" "acc" > Bact_120.tsv

fi

'''

}

process FinalOutput {

label 'big_mem'
label 'small_cpus'

input:
path "Result.fasta"
path "Result.tsv"
path "Bact_120.tsv"

output:
path "plaSquid_result.tsv"
path "plaSquid_result.fasta"

publishDir "$params.outdir/", mode: "copy"

script:

"""

Final_output.R Result.tsv Result.fasta Bact_120.tsv

"""

}


process DownPLSDB {

label 'big_mem'
label 'big_cpus'

  output:
  path "plsdb.fna"
  
  script:
  """
  
  wget https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2
  bzip2 -d plsdb.fna.bz2
  
    
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
