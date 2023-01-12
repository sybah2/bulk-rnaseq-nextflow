#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the channel for single and paired end reads
if (params.isPaired) {
    samples_ch = Channel.fromFilePairs("${params.reads}/*_{1,2}*")
}
else {
    samples_ch = Channel.fromPath("${params.reads}/*")
}
samples_ch.view()

// FastQC quality control process defination
process QualityControl{
    //container = 'biocontainers/fastqc:v0.11.9_cv7'

    publishDir "$params.results/QC_reports", mode: 'copy'

    input:
    tuple val(sampleName), path(reads)

    output:
    path "${sampleName}_QC_report"

    script:
    //sample_id = reads.getSimpleName()

    """
    FastQC1.sh "$reads" "${sampleName}"
    """
}


// Paired end trimming (with Trimmomatic) process defination 
process PaireEndTrimming {
    container = 'veupathdb/shortreadaligner'

    //this is not working - not sure why
    //publishDir "$params.results", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq"), path("${sample_id}_paired_2.fq"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    
    """ 
    
    PairedEndTrimming.sh ${projectDir} ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq  ${sample_id}_unpaired_1.fq ${sample_id}_paired_2.fq ${sample_id}_unpaired_2.fq ${params.adapters} ${sample_id}_Trimlog.txt
   
    """
}

// I haven't modified this process yet
// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {
    container = 'veupathdb/shortreadaligner'

    tag {sample_id}
    
    publishDir "${params.results}/Trimmed_reads/${sample_id}", mode: 'copy'

    input:
    tuple path(reads), val(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    """
    SingleEndTrimming.sh ${projectDir} ${reads} ${sample_id}_trim.fq ${adapters} ${sample_id}_Trimlog.txt
    
    """
}

process HisatMapping{
    container = 'veupathdb/shortreadaligner'

    tag {sample_id}

    //publishDir "$params.results/${sample_id}", mode='copy'

    input:
    tuple val(sample_id), path("${sample_id}_paired_1.fq"), path("${sample_id}_paired_2.fq")


    output:
    path("${sample_id}_sorted.bam")

    script:

    """
    hisat2 -x params.index --max-intronlen params.intronLenght -1 "${sample_id}_paired_1.fq" -2 "${sample_id}_paired_2.fq"i 2>hisat2.log  | samtools view - bS -  | samtools sort - > ${sample_id}_sorted.bam
    """

}

// work flow difination
workflow {
   // call the QC step
    fastqc = QualityControl(samples_ch)
   
   // Trimming based on paired or not
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(samples_ch)
        trimm.trimmed_fastqs.view()
        //trimm.trimm_log.view()
    }
//    else 
//    {
//        trimm = SingleEndTrimming(singleEndReads)
//        trimm.trimmed_fastqs.view()
//        //trimm.trimm_log.view()

//    }

    hisat = HisatMapping(trimm.trimmed_fastqs)
    hisat.view()
}
