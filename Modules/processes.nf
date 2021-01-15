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
    
    Rscript $baseDir/bin/Parse_paf.R ${contig_id}.paf ${contig_id}.fst

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