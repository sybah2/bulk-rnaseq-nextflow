#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { rna_seq } from  './modules/bulkRnaSeq.nf'


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------

if (params.isPaired) {
    isPairedEnd = 1
} else {
    isPairedEnd = 0
}

if (params.isStranded) {
    strandSpecific = 1
}else {
    strandSpecific = 0
}


 if(!params.reads) {
    throw new Exception("Missing parameter params.reads")
  }
if(!params.reference) {
    throw new Exception("Missing parameter params.reference")
  }

if(!params.annotation) {
    throw new Exception("Missing parameter params.annotation")
  }

if(!params.annotation) {
    throw new Exception("Missing parameter params.intronLenght")
  }


    reads_qc = Channel.fromPath("${params.reads}/*", checkIfExists: true) 

// read_ch for trimming and mapping generated based on weather the sequencing is single or paired end.
    if (params.isPaired){
        reads_ch = Channel.fromFilePairs([params.reads + '/*_{1,2}.fastq', params.reads + '/*_{1,2}.fastq.gz', params.reads + '/*_{1,2}.fq.gz'])
    } else {
        reads_ch = Channel.fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])

}

//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
    rna_seq(reads_qc, reads_ch)
}