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

// Download paired reads
process sraDownloadPaired {
    container = 'veupathdb/bowtiemapping'
    input:
    val(sra)

    output:
    tuple val(sra), path("${sra}*"), emit: trim_sample
    path("${sra}*"), emit: qc_sample

    script:
    template 'sraPairedDownload.bash'
}
// Download single reads
process sraDownloadSingle {
    container = 'veupathdb/bowtiemapping'
    input:
    val(sra)

    output:
    path("${sra}*")

    script:

    template 'sraSingleDownload.bash'
    
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

    container = 'veupathdb/shortreadaligner'

    input:
    tuple val(sample_id), path(fastqc_out)

    output:
    path("quality_check_out")

    script:

    template 'fastqc_check.bash'

}
 
 // Paired end trimming process
process paireEndTrimming {

   //container = 'veupathdb/shortreadaligner'

    input:
    path(quality_check_out)
    tuple val(sample_id), path(reads)


    output:
    val(sample_id), emit: sampleID
    path("${sample_id}_paired*.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    template 'trimmingPaired.bash'
}
// Single end process
process singleEndTrimming {

   // container = 'veupathdb/shortreadaligner'

    input:
    path(quality_check_out)
    path(reads)
    

    output:
    val(sample_id), emit: sampleID
    path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    template 'trimmingSingle.bash'

}
// Hisat mapping process for paired end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
process hisatMappingPairedEnd{
    container = 'veupathdb/shortreadaligner'

    input:
    path(quality_check_out)
    tuple path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'

    output:
    path("${sample_id}.bam")

    script:
    sample_id = paired1.getBaseName()
    if (params.intronLenght < 20) {
        template 'hisat2PairedNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hisat2Paired.bash'
    }

}

// Hisat mapping process for single end reads, taking into account weather the need for splice aware mapping or not, output coordiate sorted bam file
process hisatMappingSingleEnd{
    container = 'veupathdb/shortreadaligner'

    input:
    path(quality_check_out)
    path(paired1)
    val index 
    path 'genomeIndex.*.ht2'

    output:
    path("${sample_id}.bam")

    script:
    sample_id = paired1.getBaseName()

    if (params.intronLenght < 20) {
        template 'hisat2SingleNoSplicing.bash'
    } else if (params.intronLenght >= 20) {
        template 'hisat2Single.bash'
    }
    
}
process sortSam {
    container = 'veupathdb/shortreadaligner'
    input:
    path(sam)

    output:
    tuple val("${sample_base}"), path("*bam")

    script:
    split_name = sam.getBaseName()
    sample_base = sam.getSimpleName()

    template 'samSorting.bash'
    
}

process mergeSams {
    container = 'veupathdb/shortreadaligner'
    input:
    tuple val(sampleID), path("*.bam")
    

    output:
    path("${sampleID}.bam"), emit: sam
    tuple val(sampleID), path("*bam"), emit: bam

    script:

    template 'samMerge.bash'
    
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

    container = 'veupathdb/shortreadaligner'

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
    container = 'veupathdb/shortreadaligner'

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
        index_ch = hisatIndex(params.reference)
 
        if (params.local && params.isPaired)
        {   
            fastqc = qualityControl(reads_qc)

            chech_fastq = fastqcCheck(fastqc) | first()

            trim = paireEndTrimming(chech_fastq,reads_ch)

            reads = trim.trimmed_fastqs
                    .splitFastq( by : params.splitChunk, pe: true, file:true  )

            hisat = hisatMappingPairedEnd(chech_fastq,reads,  index_ch.genome_index_name, index_ch.ht2_files)

        } else if(!params.local && params.isPaired) 
        {
            sample = sraDownloadPaired(reads_ch) 

            sample_qc = sample.qc_sample | flatten()
            
            fastqc = qualityControl(sample_qc)

            chech_fastq = fastqcCheck(fastqc) | first()

            
            trim = paireEndTrimming(chech_fastq, sample.trim_sample)

            reads = trim.trimmed_fastqs
                    .splitFastq( by : params.splitChunk, pe: true, file:true  )

            hisat = hisatMappingPairedEnd(chech_fastq,reads,  index_ch.genome_index_name, index_ch.ht2_files)

    
        } else if(!params.local && !params.isPaired) {
            
            sample = sraDownloadSingle(reads_ch)
/*
            fastqc = qualityControl(sample)

            chech_fastq = fastqcCheck(fastqc) | first()

            trim = singleEndTrimming(chech_fastq,sample)

            reads = trim.trimmed_fastqs
                    .splitFastq( by : params.splitChunk, file:true  )

            hisat = hisatMappingSingleEnd(chech_fastq,reads, index_ch.genome_index_name, index_ch.ht2_files)
*/
        }
        else 
        {
            fastqc = qualityControl(reads_ch)

            chech_fastq = fastqcCheck(fastqc) | first()

            trim = singleEndTrimming(chech_fastq,reads_ch)
            
            reads = trim.trimmed_fastqs
                    .splitFastq( by : params.splitChunk, file:true  )

            hisat = hisatMappingSingleEnd(chech_fastq,reads, index_ch.genome_index_name, index_ch.ht2_files)
    
    }
/*
        sortedsam = sortSam(hisat) 

        samSet = sortedsam.groupTuple(sort: true)

        mergeSam = mergeSams(samSet)

        sortedByName = sortBams(mergeSam.bam)

        hisatCount = htseqCounting(sortedByName, params.annotation)

        beds_stats = bedBamStats(mergeSam.bam)

        spliceCounts = spliceCrossingReads(mergeSam.bam)
*/
}
