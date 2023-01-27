#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// The processes defined below using template bash files located in the templates folder in workflow directory. 

// Define the channel for single and paired end reads, the second channel for the paired is just for the QC step
// read_qc used just for the QC step
reads_qc = Channel.fromPath("${params.reads}/*", checkIfExists: true) 

// read_ch for trimming and mapping generated based on weather the sequencing is single or paired end.
if (params.isPaired){
    reads_ch = Channel.fromFilePairs([params.reads + '/*_{1,2}.fastq', params.reads + '/*_{1,2}.fastq.gz', params.reads + '/*_{1,2}.fq.gz'])
} else {
    reads_ch = Channel.fromPath([params.reads + '/*.fastq', params.reads + '/*.fastq.gz', params.reads + '/*.fq.gz'])

}

// Hist indexing process. 
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

// FastQC quality control process
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

// Fastqc_check process to get the phred encoding
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
 
 // Paired end trimming process
process PaireEndTrimming {
    tag {sample_id}

    container = 'saikou'

    input:
    path("quality_check_out")
    tuple val(sample_id), path(reads)


    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    template 'trimming_paired.bash'
}
// Single end process
process SingleEndTrimming {
    tag {sample_id}

    container = 'saikou'

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
// Hisat mapping process for paired end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
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
    if (params.intronLenght < 20) {
        template 'hista2PairedNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hista2Paired.bash'
    }

}

// Hisat mapping process for single end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
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
    if (params.intronLenght < 20) {
        template 'hista2SingleNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hista2Single.bash'
    }
    

}

// Process to sort bam file by names
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

// HTSeq counting process, takes into account weather the sequencing is stranded or not. 
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
// Process to generate splice juctions
process SpliceCrossingReads{   

    container = 'saikou' 
    tag {sample_id}

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("junctions.tab")

    script:
    template 'spliceCorssReads.bash'
}

// Generate bam statistic and bed files from the bam files.
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
    fastqc = QualityControl(reads_qc)
    chech_fastq = Fastqc_check(fastqc)
    
    index_ch = HisatIndex(params.reference)
    
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(chech_fastq,reads_ch)
        hisat = HisatMappingPairedEnd(chech_fastq,trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
    }
    else 
    {
        trimm = SingleEndTrimming(chech_fastq,reads_ch)
        hisat = HisatMappingSingleEnd(chech_fastq,trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
    
    }
 
    sortedByName = SortBams(hisat)
    hisatCount = HtseqCounting(sortedByName)
    beds_stats = BedBamStats(hisat)
    spliceCounts = SpliceCrossingReads(hisat)

}
