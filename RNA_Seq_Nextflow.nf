#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the channel for single and paired end reads, the second channel for the paired is just for the QC step
if (params.isPaired){
read_pairs_ch = Channel.fromFilePairs("${params.reads}/*_{1,2}*", checkIfExists: true)
singleEndReads = Channel.fromPath("${params.reads}/*", checkIfExists: true)
//singleEndReads.view()
} else {
    singleEndReads = Channel.fromPath("${params.reads}/*", checkIfExists: true)
}

process HisatIndex {

    tag 'Hisat2 indexing'

    container = 'veupathdb/shortreadaligner'

    publishDir "${projectDir}/index/", mode: 'copy'    

    input:
    path(reference)

    output: 
    path 'genomeIndex*.ht2', emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

    script:
        template 'hista2_index.bash' 
}

// FastQC quality control process defination
process QualityControl {

    tag {sample_id}

    container = 'biocontainers/fastqc:v0.11.9_cv7'
    input:
    path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_QC")

    script:
    sample_id = reads.getSimpleName()
    template 'fastqc.bash'
}

process Fastqc_check {

    tag {sample_id}

    container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sample_id), path(fastqc_out)

    output:
    path("quality_check_out")

    script:

    template 'fastqc_check.bash'

}
// Paired end trimming (with Trimmomatic) process defination 
process PaireEndTrimming {
    tag {sample_id}

    input:
    path("quality_check_out")
    tuple val(sample_id), path(reads)
    
    

    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    template 'trimming_paired.bash'
}

// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {
    tag {sample_id}

    input:
    path("quality_check_out")
    path(reads)
    

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    template 'trimming_single.bash'

}

process HisatMappingPairedEnd{
    container = 'veupathdb/shortreadaligner'

    tag {sample_id}

    input:
    path("quality_check_out")
    tuple val(sample_id), path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    
    template 'hista2Paired.bash'

}

process HisatMappingSingleEnd{
    tag {sample_id}
    container = 'veupathdb/shortreadaligner'

    input:
    path("quality_check_out")
    tuple val(sample_id), path(paired1)
    val index 
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    template 'hista2Single.bash'

}

process SortBams{
    tag {sample_id}
    container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sortedByName.bam")

    script:
    template 'BamSortbyName.bash'

}
process HtseqCounting{
    tag {sample_id} 
    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'genomicpariscentre/htseq'
      
    input:
    tuple val(sample_id), path(bam)

    output:
    path("*.counts")
    
    script:
    if (params.isStranded) {
        template 'htseqCounStranded.bash'
    } else {
        template 'htseqCounUnStranded.bash'
    }
    
}

process SpliceCrossingReads{   
    tag {sample_id}

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'saikou' 
    input:
    tuple val(sample_id), path(bam)

    output:
    path("junctions.tab")

    script:
    template 'spliceCorssReads.bash'
}

process BedBamStats{
    tag {sample_id}

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'saikou'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("*.bed"), path("mappingStats.txt")
    
    script:

    if (params.isPaired && params.isStranded){
        template 'BedFileStatsStrandedPaired.bash'
    } else if  (!params.isPaired && params.isStranded) {
        template 'BedFileStatsStrandedUnpaired.bash'
    } else if  (params.isPaired && !params.isStranded) {
        template 'BedFileStatsUnstrandedPaired.bash'
    } else if  (!params.isPaired && !params.isStranded) {
        template 'BedFileStatsUnstrandedUnpaired.bash'
    }
    
}
// work flow difination
workflow {
   // call the QC step
    fastqc = QualityControl(singleEndReads)
    // fastqc check
    chech_fastq = Fastqc_check(fastqc)
    
    index_ch = HisatIndex(params.reference)
     
   // Trimming based on paired or not
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(chech_fastq,read_pairs_ch)
        hisat = HisatMappingPairedEnd(chech_fastq,trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
        //hisat.view()
    }
    else 
    {
        trimm = SingleEndTrimming(chech_fastq,singleEndReads)
        hisat = HisatMappingSingleEnd(chech_fastq,trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
        //hisat.view()
    }
 
    sortedByName = SortBams(hisat)
    hisatCount = HtseqCounting(sortedByName)
    beds_stats = BedBamStats(hisat)
    spliceCounts = SpliceCrossingReads(hisat)

}
