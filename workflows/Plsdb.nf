#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.plsdbURL = "https://ndownloader.figshare.com/files/23582252"
params.mmi = "$baseDir/plsdb.mmi"

include {DownPLSDB} from '../Modules/processes.nf'
include {FormtPLSDB} from '../Modules/processes.nf'

workflow SetPlsdb {

main:

plsdbmmi = file( params.mmi )

if ( plsdbmmi.exists() ) {

Channel.value( plsdbmmi )
       .set{ FmtPlsdb_ch }

} else {

         DownPLSDB()
         DownPLSDB.out
                  .set{ plsdb_ch }
      
        FormtPLSDB( plsdb_ch )
        FormtPLSDB.out
                  .set{ FmtPlsdb_ch }

}


emit:
dbs_ch = FmtPlsdb_ch

}


