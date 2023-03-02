#!/usr/bin/env nextflow
nextflow.enable.dsl=2
import nextflow.splitter.CsvSplitter

def fetchRunAccessions( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> run_accessions = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       run_accessions.add( row['run_accession'] )
    }
    return run_accessions
}

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { rna_seq } from  './modules/bulkRnaSeq.nf'


//---------------------------------------------------------------
// Param Checking 
//---------------------------------------------------------------


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

// read_ch for trimming and mapping generated based on weather the sequencing is single or paired end. If the input is SRA sample IDs they should be saved in csv file. 
    if (params.local && params.isPaired){
        reads_ch = Channel.fromFilePairs([params.reads + '/*_{1,2}.fastq', params.reads + '/*_{1,2}.fastq.gz', params.reads + '/*_{1,2}.fq.gz'])
    } else if (!params.local) {
        input = fetchRunAccessions(params.sraAccession)
        reads_ch = Channel.fromList(input)
    }
    else {
        reads_ch = Channel.fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])

}

//--------------------------------------
// Process the workflow
//-------------------------------------

workflow {
    rna_seq(reads_qc, reads_ch)
}