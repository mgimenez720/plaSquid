#!/usr/bin/env nextflow


nextflow.preview.dsl = 2


include { Renamecntgs } from "../Modules/processes.nf"
include { AnnotContigs } from "../Modules/processes.nf"
include { MobSearch } from "../Modules/processes.nf"
include { MobFilter } from "../Modules/processes.nf"
include { GeneRetrieve } from "../Modules/processes.nf"
include { MobExtract } from "../Modules/processes.nf"


workflow MOBextract {
take:

fasta_ch

main:

Renamecntgs(fasta_ch)

Renamecntgs.out
           .set{ assembl_ch }

AnnotContigs( assembl_ch )

AnnotContigs.out
            .set{ search_ch }

//MOB search

MobSearch( search_ch )

MobSearch.out
         .set{ Mob_ch }
MobFilter( Mob_ch )

MobFilter.out
         .set{ mob_ret_ch }

MobExtract( search_ch, mob_ret_ch )
MobExtract.out
          .set{ mob_extr_ch }

emit:
mob_extr_ch

}






