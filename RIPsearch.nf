#!/usr/bin/env nextflow

nextflow.preview.dsl = 2

params.contigs = "*.fasta"

include { AnnotContigs } from "./Modules/processes.nf"
include { RepSearch } from "./Modules/processes.nf"
include { DomainArch } from "./Modules/processes.nf"
include { IncSearch } from "./Modules/processes.nf"
include { RnaSearch } from "./Modules/processes.nf"
include { IncClassif } from "./Modules/processes.nf"
include { IncFilter } from "./Modules/processes.nf"
include { MobSearch } from "./Modules/processes.nf"
include { MobFilter } from "./Modules/processes.nf"
include { GeneRetrieve } from "./Modules/processes.nf"


workflow {

Channel
  .fromPath(params.contigs, checkIfExists: true)
  .ifEmpty { exit 1, "Non fasta files found: ${params.contigs}" }
  .set{ fasta_ch }

AnnotContigs(fasta_ch)

AnnotContigs.out
            .set{ search_ch }

IncSearch(search_ch)
MobSearch(search_ch)

IncSearch.out
         .set{ Inc_ch }

RnaSearch(search_ch)

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

GeneRetrieve(Inc_ret_ch, Mob_ret_ch, search_ch)

//Continuar trabajando con los dominios para lograr un modo comprehensive


}
