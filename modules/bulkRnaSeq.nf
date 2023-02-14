#!/usr/bin/env nextflow
nextflow.enable.dsl=2


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

// Hisat indexing process. 
process hisatIndex {

    container = 'veupathdb/shortreadaligner'

    input:
    path(reference)

    output: 
    path 'genomeIndex*.ht2', emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

    script:
        if (params.createIndex) {
            template 'hisat2Index.bash' 
        } else {
            template 'hisat2NoIndex.bash'
        }
    stub:
    """
    touch genomeIndex.1.ht2
    """
}

// fastQC quality control process

process qualityControl {
    container = 'biocontainers/fastqc:v0.11.9_cv7'
    
    input:
    path reads

    output:
    tuple val(sample_id), path("${sample_id}_fastqc_out")

    script:
    sample_id = reads.getSimpleName()
    template 'fastqc.bash'
}


// Fastqc_check process to get the phred encoding
process fastqcCheck {

    container = 'veupathdb/shortreadaligner:branch-saikou'

    input:
    tuple val(sample_id), path(fastqc_out)

    output:
    path("quality_check_out")

    script:

    template 'fastqc_check.bash'

}
 
 // Paired end trimming process
process paireEndTrimming {

   container = 'veupathdb/shortreadaligner:branch-saikou'

    input:
    path(quality_check_out)
    tuple val(sample_id), path(reads)


    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    template 'trimmingPaired.bash'
}
// Single end process
process singleEndTrimming {

    container = 'veupathdb/shortreadaligner:branch-saikou'

    input:
    path(quality_check_out)
    path(reads)
    

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    template 'trimmingSingle.bash'

}
// Hisat mapping process for paired end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
process hisatMappingPairedEnd{
    container = 'veupathdb/shortreadaligner:branch-saikou'

    input:
    path(quality_check_out)
    tuple val(sample_id), path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    if (params.intronLenght < 20) {
        template 'hisat2PairedNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hisat2Paired.bash'
    }

}

// Hisat mapping process for single end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
process hisatMappingSingleEnd{
    container = 'veupathdb/shortreadaligner:branch-saikou'

    input:
    path(quality_check_out)
    tuple val(sample_id), path(paired1)
    val index 
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    if (params.intronLenght < 20) {
        template 'hisat2SingleNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hisat2Single.bash'
    }
    
}

// Process to sort bam file by names
process sortBams{
    container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sortedByName.bam")

    script:
    template 'BamSortbyName.bash'

}

// HTSeq counting process, takes into account weather the sequencing is stranded or not. 
process htseqCounting{
    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'
      
    input:
    tuple val(sample_id), path(bam)
    path(annotation)

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
process spliceCrossingReads{   

    container = 'veupathdb/shortreadaligner:branch-saikou'

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("junctions.tab")

    script:
    template 'spliceCrossReads.bash'
}

// Generate bam statistic and bed files from the bam files.
process bedBamStats{
    container = 'veupathdb/shortreadaligner:branch-saikou'

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("*.bed"), path("mappingStats.txt")
    
    script:
        template 'BedFileStats.bash'
}
// work flow difination
workflow rna_seq {

    take:
        reads_qc
        reads_ch

    main: 
        fastqc = qualityControl(reads_qc)
        chech_fastq = fastqcCheck(fastqc)
  
        index_ch = hisatIndex(params.reference)
 
        if (params.isPaired) 
        {
            trim = paireEndTrimming(chech_fastq,reads_ch)
            hisat = hisatMappingPairedEnd(chech_fastq,trim.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
        }
        else 
        {
            trim = singleEndTrimming(chech_fastq,reads_ch)
            hisat = hisatMappingSingleEnd(chech_fastq,trim.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
    
    }

        sortedByName = sortBams(hisat)
        hisatCount = htseqCounting(sortedByName, params.annotation)
        beds_stats = bedBamStats(hisat)
        spliceCounts = spliceCrossingReads(hisat)

}
