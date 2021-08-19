#!/usr/bin/env nextflow

nextflow.preview.dsl = 2


include { Renamecntgs } from "../Modules/processes.nf"
include { AnnotContigs } from "../Modules/processes.nf"
include { RepSearch } from "../Modules/processes.nf"
include { DomainArch } from "../Modules/processes.nf"
include { FilterDom } from "../Modules/processes.nf"
include { IncSearch } from "../Modules/processes.nf"
include { RnaSearch } from "../Modules/processes.nf"
include { IncClassif } from "../Modules/processes.nf"
include { IncFilter } from "../Modules/processes.nf"
include { GeneRetrieve } from "../Modules/processes.nf"
include { RipExtract } from "../Modules/processes.nf"

workflow RIPextract {
take:

fasta_ch

main:

Renamecntgs(fasta_ch)

Renamecntgs.out
           .set{ assembl_ch }

AnnotContigs( assembl_ch )

AnnotContigs.out
            .set{ search_ch }


//Domain architecture RIP search
RepSearch( search_ch )

RepSearch.out
         .set{ domar_ch }

DomainArch( domar_ch )

DomainArch.out
          .set{ domfil_ch }

FilterDom( domfil_ch )

FilterDom.out
         .set{ repdom_ch }


//Inc group classification
IncSearch( search_ch )

IncSearch.out
         .set{ Inc_ch }

RnaSearch( assembl_ch )

RnaSearch.out
         .set{ Rna_ch }

IncClassif(Inc_ch, Rna_ch)

IncClassif.out
          .set{ Filter_ch }

IncFilter( Filter_ch )
IncFilter.out
         .set{ Inc_ret_ch }

//Collapse information and report output

RipExtract( search_ch, Inc_ret_ch, repdom_ch)
RipExtract.out
          .set{ rip_extr_ch }

emit:
rip_extr_ch

}
