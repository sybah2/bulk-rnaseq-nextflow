#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads = "$projectDir/Fastq"
params.index = "$projectDir/index/pfal3D7"
params.isPaired = false
params.isStranded = true
params.intronLenght = 3636
params.adapters = "$projectDir/All_adaptors-PE.fa"
params.results = 'Results'

// Define the channel for single and paired end reads
read_pairs_ch = Channel.fromFilePairs("${params.reads}/*_{1,2}*")
singleEndReads = Channel.fromPath("${params.reads}/*")

// FastQC quality control process defination
process QualityControl{

    tag {sample_id}

    publishDir "${params.results}/QC_reports/", mode: 'copy'


    input:
    path(reads)

    output:
    path "${sample_id}_QC_report"

    script:
    sample_id = reads.getSimpleName()

    """
    FastQC1.sh "$reads", "${sample_id}"
    """
}


// Paired end trimming (with Trimmomatic) process defination 
process PaireEndTrimming {

    tag {sample_id}

    publishDir "${params.results}/Trimmed_reads/${sample_id}", mode: 'copy'

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

// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {

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
    touch ${sample_id}_trim1.fq
    SingleEndTrimming.sh ${projectDir} ${reads} ${sample_id}_trim.fq ${params.adapters} ${sample_id}_Trimlog.txt
    
    """

}

process HisatMapping{
    
    container = 'veupathdb/shortreadaligner'


    publishDir "${params.results}/Bam/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1), path(paired2)
    


    output:
    path("${sample_id}_sorted.bam")

    script:
    
    """
    hisat2 -x params.index --max-intronlen params.intronLenght 20 ${paired1} -2 ${paired2} 2>hisat2.log  | samtools view - bS -  | samtools sort - > ${sample_id}_sorted.bam
    """

}

// work flow difination
workflow {
   // call the QC step
    fastqc = QualityControl(singleEndReads)
   
   // Trimming based on paired or not
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(read_pairs_ch)
        trimm.trimmed_fastqs.view()
        //trimm.trimm_log.view()
    }
    else 
    {
        trimm = SingleEndTrimming(singleEndReads)
        trimm.trimmed_fastqs.view()
        //trimm.trimm_log.view()

    }

    hisat = HisatMapping(trimm.trimmed_fastqs)
    //hisat.view()
}
