#!/usr/bin/env nextflow

nextflow.preview.dsl = 2


include { Splitter } from "../Modules/processes.nf"
include { Mapping_pr } from "../Modules/processes.nf"
include { Parse_paf } from "../Modules/processes.nf"
include { RetrievePlasmids } from "../Modules/processes.nf"

workflow Minidist {

take:
fasta_ch
dbs_ch

main:

Splitter(fasta_ch)
   
   Splitter.out
           .set{ fasta_spl_ch }

Mapping_pr(fasta_spl_ch, dbs_ch)

Mapping_pr.out 
          .set{ paf_parse_ch }

Parse_paf(paf_parse_ch, fasta_ch)   
Parse_paf.out
         .set{minidist_ch}



emit:
minidist_ch

}



