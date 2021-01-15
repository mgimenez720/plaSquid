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


process MakeBldb {

   input:
   path contigs

   output:
   path "blast.db*"

   script:
   """
   
   awk '/^>/{print ">Contig-" ++i; next}{print}' ${contigs} > tmp

   makeblastdb -in tmp -dbtype nucl -out blast.db -parse_seqids -hash_index
   
   
   """
}

process RetrievePlasmids {

   input:
   path "*.tsv"
   path contigs
   

   output:
   path "contigs_list.txt"
   path "plasmids.fasta"

   script:
   """
   cat *.tsv | cut -f1 > contigs_list.txt

   awk '/^>/{print ">Contig-" ++i; next}{print}' ${contigs} > tmp

   makeblastdb -in tmp -dbtype nucl -out blast.db -parse_seqids -hash_index
   
   blastdbcmd -db blast.db -dbtype nucl -entry_batch contigs_list.txt -out plasmids.fasta

   """


}