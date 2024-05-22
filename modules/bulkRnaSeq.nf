#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


if (params.isStranded) {
    strandSpecific = 1
}else {
    strandSpecific = 0
}
if (params.isPaired) {
    isPairedEnd = 1
} else {
    isPairedEnd = 0
}

/*
This procees down the fastq file from the sequence read archive if needed

@val is the accession number for he fastq file

Output is the fastq file
*/

process downloadFiles {
  container = 'veupathdb/bowtiemapping'
  
  input:
    val id
    
  output:
    tuple val(id), path("${id}**.fastq"), emit: samples

  script:
    template 'downloadFiles.bash'
}

/*
This pocess create the index file of the reference genome using hisat2

*/

process createIndex {
  container = 'veupathdb/shortreadaligner'

  input:
    val(organismAbbv)
    path(reference)

  output: 
    path "${organismAbbv}*.ht2", emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

  script:
    template 'createIndex.bash'
}

/*
If the index is already generated this process copy the index files to the right location

*/
process copyIndex {
  container = 'veupathdb/shortreadaligner'

  input:
    val(organismAbbv)
    path(hisat2Index)

  output: 
    path "${organismAbbv}*.ht2", emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

  script:
    template 'copyIndex.bash'
}

/*
Perform the quality control of the fastq files
*/

process qualityControl {
  container = 'biocontainers/fastqc:v0.11.9_cv7'
    
  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_fastqc_out"), emit: qcOutput
    tuple val(sample_id), path(reads), emit: sample

  script:
    template 'fastqc.bash'
}

/*
This process does the fastqc check to determine the QC for each read to be used for the trimming
*/

process fastqcCheck {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sample_id), path(fastqc_out)

  output:
    path("quality_check_out")

  script:
    template 'fastqc_check.bash'
}

/*
Quality trimming for paired end reads
*/
process pairedEndTrimming {
  container = 'veupathdb/shortreadaligner'

  input:
    path(quality_check_out)
    tuple val(sample_id), path(reads)

  output:
    val(sample_id), emit: sampleID
    path("${sample_id}.{1,2}.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
  script:
    template 'trimmingPaired.bash'
}

/*
Quality trimming for single end reads
*/
process singleEndTrimming {
  container = 'veupathdb/shortreadaligner'

  input:
    path(quality_check_out)
    tuple val(sample_id), path(reads)

  output:
    val(sample_id), emit: sampleID
    path("${sample_id}.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

  script:
    sample_id = reads.getSimpleName()
    template 'trimmingSingle.bash'
}

/*
Hisat mapping for paired end reads. This process split each fastq files for faster mapping
*/
process hisatMappingPairedEnd{
  container = 'veupathdb/shortreadaligner'

  input:
    path(quality_check_out)
    tuple path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'
    val intronLength

  output:
    path("${sample_id}.bam")

  script:
    sample_id = paired1.getBaseName()
    template 'hisatMappingPairedEnd.bash'
}

/*
Hisat mapping for single end reads. This process split each fastq files for faster mapping
*/

process hisatMappingSingleEnd{
  container = 'veupathdb/shortreadaligner'

  input:
    path(quality_check_out)
    path(read)
    val index 
    path 'genomeIndex.*.ht2'
    val intronLength

  output:
    path("${sample_id}.bam")

  script:
    sample_id = read.getBaseName()
    template 'hisatMappingSingleEnd.bash'
}

/*
This process sort the alignment file
*/

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

/*
This process merge the alignment file for each sample into one bam file
*/
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

/*
This process sort the merge bam file
*/
process sortBams{
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path("${sample_id}_sortedByName.bam")

  script:
    template 'BamSortbyName.bash'
}
/*
This process does the counting of the reads that map onto each genes
*/

process htseqCounting{
  publishDir "${params.results}/${sample_id}", mode: 'copy'

  container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'
      
  input:
    tuple val(sample_id), path(bam)
    path(annotation)
    val(isCds)
    val isStranded

  output:
    path("*.counts")
    
  script:
    template 'htseqCounting.bash'    
}
/*
This process determine splice crossing reads
*/

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

/*
The process generate the statistics for the alignment file
*/
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


workflow rna_seq {

  take:
    reads_ch

  main:

      if (params.createIndex){
        index_ch = createIndex(params.organismAbbv, params.reference)
      } else{
        index_ch = copyIndex(params.organismAbbv, params.hisat2Index)
      }


    if (params.local) {
      fastqc = qualityControl(reads_ch)
    }
    else {
      sample = downloadFiles(reads_ch)
      fastqc = qualityControl(sample.samples)
    }

    check_fastq = fastqcCheck(fastqc.qcOutput) | first()

    if (params.isPaired) {
        trim = pairedEndTrimming(check_fastq,fastqc.sample)
        reads = trim.trimmed_fastqs.splitFastq( by : params.splitChunk, pe: true, file:true)
        hisat = hisatMappingPairedEnd(check_fastq,reads, index_ch.genome_index_name, index_ch.ht2_files, params.intronLength) 
    }
    else {
        trim = singleEndTrimming(check_fastq,fastqc.sample)
        reads = trim.trimmed_fastqs.splitFastq( by : params.splitChunk, file:true  )
        hisat = hisatMappingSingleEnd(check_fastq,reads, index_ch.genome_index_name, index_ch.ht2_files, params.intronLength)
    }

   sortedsam = sortSam(hisat) 
   samSet = sortedsam.groupTuple(sort: true)
   mergeSam = mergeSams(samSet)
   sortedByName = sortBams(mergeSam.bam)
   hisatCount = htseqCounting(sortedByName, params.annotation, params.isCds, params.isStranded)
   beds_stats = bedBamStats(mergeSam.bam)
   spliceCounts = spliceCrossingReads(mergeSam.bam)

}