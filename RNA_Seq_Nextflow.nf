#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the channel for single and paired end reads
read_pairs_ch = Channel.fromFilePairs("${params.reads}/*_{1,2}*", checkIfExists: true)
singleEndReads = Channel.fromPath("${params.reads}/*", checkIfExists: true)
//singleEndReads.view()

process HisatIndex {
    container = 'veupathdb/shortreadaligner'

    publishDir "${projectDir}/index/", mode: 'copy'    

    input:
    path(reference)

    output: 
    path 'genomeIndex*.ht2', emit: ht2_files
    val 'genomeIndex' , emit: genome_index_name

    script:

    """
    hisat2-build ${reference} genomeIndex
    """  
}

// FastQC quality control process defination
process QualityControl {
    publishDir "${params.results}/QC_reports/", mode: 'copy'

    input:
    path(singleEndReads)

    output:
    // tuple file('*.html'), file('*.zip') 
    path("${sample_id}_QC")


    script:
    sample_id = singleEndReads.getSimpleName()

    """
    mkdir ${sample_id}_QC
    fastqc ${params.fastqcExtension}  -o "${sample_id}_QC" --extract
    """
}

// Paired end trimming (with Trimmomatic) process defination 
process PaireEndTrimming {

    tag {sample_id}

    publishDir "${params.results}/Trimmed_reads/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    """  
    PairedEndTrimming.sh ${projectDir} ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz ${params.adapters} ${sample_id}_Trimlog.txt 
    """
}

// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {

    tag {sample_id}
    
    publishDir "${params.results}/Trimmed_reads/${sample_id}", mode: 'copy'

    input:
    tuple path(reads), val(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    """
    touch ${sample_id}_trim.fq.gz
    SingleEndTrimming.sh ${projectDir} ${reads} ${sample_id}_trim.fq.gz ${params.adapters} ${sample_id}_Trimlog.txt  
    """

}

process HisatMappingPairedEnd{

    publishDir "${params.results}/Bam/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1), path(paired2)
    //val hisat2_index 
    //path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    
    """
    Hisat2PairedMapping.sh ${params.index} ${params.intronLenght} ${paired1} ${paired2} ${sample_id}
    """

}

process HisatMappingSingleEnd{

    publishDir "${params.results}/Bam/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1)
    //val hisat2_index 
    //path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    """
    Hisat2SingleEndMapping.sh ${params.index} ${params.intronLenght} ${paired1} ${sample_id}  
    """

}

process SortBams{

    publishDir "${params.results}/Bam/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sortedByName.bam")

    """
    samtools sort -n ${bam} -o "${sample_id}_sortedByName.bam"
    """

}

process HtseqCountingUnstranded{
    
    publishDir "${params.results}/Counts/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("genes.htseq-union.unstranded.counts"), path("genes.htseq-union.unstranded.nonunique.counts")
    
    script:

    """
    UnstrandedCounting.sh ${bam} ${params.annotation} "genes.htseq-union.unstranded.counts" "genes.htseq-union.unstranded.nonunique.counts"

    """
}

process HtseqCountingStranded{
    
    publishDir "${params.results}/Counts/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("genes.htseq-union.firststrand.counts"), path("genes.htseq-union.secondstrand.counts"), path("genes.htseq-union.firststrand.nonunique.counts"), path("genes.htseq-union.secondstrand.nonunique.counts")
    
    script:

    """
    StrandedCounting.sh ${bam} ${params.annotation} "genes.htseq-union.firststrand.counts" "genes.htseq-union.secondstrand.counts" "genes.htseq-union.firststrand.nonunique.counts" "genes.htseq-union.secondstrand.nonunique.counts"
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

    }
    else 
    {
        trimm = SingleEndTrimming(singleEndReads)
    }
    
    index_ch = HisatIndex(params.reference)
    
    if (params.isPaired){

    hisat = HisatMappingPairedEnd(trimm.trimmed_fastqs)

    } else{

    hisat = HisatMappingSingleEnd(trimm.trimmed_fastqs)
    }
    
    sortedByName = SortBams(hisat)

    if (params.isStranded) {

        hisatCount = HtseqCountingStranded(sortedByName)

    } else
    {
        hisatCount = HtseqCountingUnstranded(sortedByName)
    }
}
