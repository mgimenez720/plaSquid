#!/bin/bash/env nextflow

//include processes
include { AnnotContigs } from "../Modules/processes.nf"
include { RepSearch } from "../Modules/processes.nf"
include { DomainArch } from "../Modules/processes.nf"
include { FilterDom } from "../Modules/processes.nf"

workflow DomSearch {

take:

fasta_ch

main:

AnnotContigs(fasta_ch)

AnnotContigs.out
            .set{ search_ch }

RepSearch( search_ch )

RepSearch.out
         .set{ domar_ch }

DomainArch( domar_ch )

DomainArch.out
          .set{ domfil_ch }

FilterDom( domfil_ch )

FilterDom.out
         .set{ repdom_ch }

emit:

repdom_ch

}