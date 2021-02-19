#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.plsdbURL = "https://ndownloader.figshare.com/files/23582252"
params.mmi = "plsdb.mmi"

include {DownPLSDB} from '../Modules/processes.nf'
include {FormtPLSDB} from '../Modules/processes.nf'

workflow SetPlsdb {

main:

plsdbmmi = file( params.mmi )

if ( ! plsdbmmi.exists() ) {

 DownPLSDB()
 DownPLSDB.out
          .set{ plsdb_ch }

} else {

 if ( plsdbmmi.getExtension() == "fasta" || "fna" ) {

   Channel.value( plsdbmmi )
          .set{ plsdb_ch }
 } else {

  println("Unrecognized plsdb extension (not .fasta nor .fna).")
  system.exit(1)

 }

} 

FormtPLSDB( plsdb_ch )
FormtPLSDB.out
        .set{ FmtPlsdb_ch }


emit:
dbs_ch = FmtPlsdb_ch

}


