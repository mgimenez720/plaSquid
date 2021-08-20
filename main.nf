#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//Parameters definition
params.contigs = "*.fasta"
params.plsdbURL = "https://ndownloader.figshare.com/files/23582252"
params.mmi = "plsdb.mmi"
params.outdir = "Results"

params.minidist = false
params.repsearch = false
params.ripextract = false

params.help = false


//Include workflows

include { SetPlsdb } from './workflows/Plsdb.nf'
include { Minidist } from './workflows/Minidist.nf'
include { RIPsearch } from './workflows/RIPsearch.nf'
include { RIPextract } from './workflows/RIPextract.nf'

//Include modules
include { SumOutput } from './Modules/processes.nf'
include { MinidistOut } from './Modules/processes.nf'
include { RepsearchOut } from './Modules/processes.nf'

def sayHi(){
  log.info '''
.           __       _____                _      __
    ____   / /____  / ___/ ____ _ __  __ (_)____/ /
   / __ ` / // __ `/ /_   / __ `// / / // // __  /
  / /_/ // // /_/ /___/ // /_/ // /_/ // // /_/ /  
 / .___//_/ '__,_//____/ `__  / '__,_//_/ ' __,_/  
/_/                        /_/                    
--------------------------------------------------
-A pipeline for Plasmids Sequences Identification-
--------------------------------------------------
'''
}

sayHi()

def helpMessage() {
    log.info """
    Usage:
    
    nextflow run plaSquid.nf --contigs 'data/*.fasta' 
    
    arguments:
    
    --contigs       Path to input data (must be surrounded with quotes).
    --mmi           Path to Minimap2 indexed (.mmi) or fasta (.fasta/.fna) plsdb database-
    --outdir        Path to output directory where results are written.

    --help          Print help message and exit 

    subworkflows:
    
    --minidist      Run mapping of contigs against plsdb database. 
    --repsearch     Run search and classification of RIP and MOB (Rel) genes.
    --ripextract    Run extraction of RIP sequences.

    profiles:
    
    -profile conda  Installs dependencies using a conda environment
    -profile server runs using 15 cpus and 50 Gb
    -profile test   tests dependencies and normal functioning


    Authors:
    
    Matías Giménez
    Ignacio Ferrés
    Gregorio Iraola


Microbial Genomics Laboratory
Institut Pasteur Montevideo (Uruguay)

 """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

workflow {

Channel
  .fromPath(params.contigs, checkIfExists: true)
  .ifEmpty { exit 1, "Non fasta files found: ${params.contigs}" }
  .set{ fasta_ch }

SetPlsdb()
SetPlsdb.out
        .set{ dbs_ch }

if (params.repsearch) {
  RIPsearch( fasta_ch )
  RIPsearch.out
           .set{ gene_search_ch }
  RepsearchOut( fasta_ch, gene_search_ch )

} else if (params.minidist) {

   Minidist( fasta_ch, dbs_ch )
   Minidist.out
        .set{ minidist_ch }
   MinidistOut( minidist_ch, fasta_ch )

} else if (params.ripextract) {

   RIPextract( fasta_ch )

} else {

  RIPsearch( fasta_ch )
  RIPsearch.out
           .set{ gene_search_ch }

  Minidist( fasta_ch, dbs_ch )
  Minidist.out
          .set{ minidist_ch }


  SumOutput( minidist_ch, gene_search_ch, fasta_ch )

}
}





