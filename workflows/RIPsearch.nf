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
include { MobSearch } from "../Modules/processes.nf"
include { MobFilter } from "../Modules/processes.nf"
include { GeneRetrieve } from "../Modules/processes.nf"


workflow RIPsearch {
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


//Inc and MOB group classification
IncSearch( search_ch )
MobSearch( search_ch )

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

MobSearch.out
         .set{ Mob_ch }
MobFilter( Mob_ch )

MobFilter.out
         .set{ Mob_ret_ch }

//Collapse information and report output

GeneRetrieve( Inc_ret_ch, Mob_ret_ch, repdom_ch, assembl_ch )
GeneRetrieve.out[1]
            .set{ gene_search_ch }


emit:
gene_search_ch

}