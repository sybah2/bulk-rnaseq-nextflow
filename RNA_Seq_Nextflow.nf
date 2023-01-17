#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the channel for single and paired end reads
read_pairs_ch = Channel.fromFilePairs("${params.reads}/*_{1,2}*", checkIfExists: true)
singleEndReads = Channel.fromPath("${params.reads}/*", checkIfExists: true)
//singleEndReads.view()

process HisatIndex {

    container = 'dceoy/hisat2'

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
    //publishDir "${params.results}/", mode: 'copy'

    input:
    path(singleEndReads)

    output:
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

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    """  
    PairedEndTrimming.sh ${projectDir} ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz ${params.adaptersPE} ${sample_id}_Trimlog.txt 
    """
}

// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {

    tag {sample_id}
    
    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    path(reads)
    val(adapters)

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    """
    touch ${sample_id}_trim.fq.gz
    SingleEndTrimming.sh ${projectDir} ${reads} ${sample_id}_trim.fq.gz ${params.adaptersSE} ${sample_id}_Trimlog.txt  
    """

}

process HisatMappingPairedEnd{

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    
    """
    hisat2 -x ${index} --max-intronlen ${params.intronLenght} -1 ${paired1} -2 ${paired2} 2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}_sorted.bam
    """

}

process HisatMappingSingleEnd{

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1)
    val index 
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}_sorted.bam")

    script:
    """
    hisat2 -x ${index} --max-intronlen ${params.intronLenght} -U ${paired1}  2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}_sorted.bam
    
    """

}

process SortBams{

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}_sortedByName.bam")

    """
    samtools sort -n ${bam} -o "${sample_id}_sortedByName.bam"
    """

}

process HtseqCountingUnstranded{
    
    //publishDir "${params.results}/${sample_id}", mode: 'copy'
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
    
    publishDir "${params.results}/${sample_id}", mode: 'copy'
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("genes.htseq-union.firststrand.counts"), path("genes.htseq-union.secondstrand.counts"), path("genes.htseq-union.firststrand.nonunique.counts"), path("genes.htseq-union.secondstrand.nonunique.counts")
    
    script:

    """
    StrandedCounting.sh ${bam} ${params.annotation} "genes.htseq-union.firststrand.counts" "genes.htseq-union.secondstrand.counts" "genes.htseq-union.firststrand.nonunique.counts" "genes.htseq-union.secondstrand.nonunique.counts"
    """
}

process SpliceCrossingReads{

    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'saikou' 
    input:
    tuple val(sample_id), path(bam)

    output:
    path("junctions.tab")

    script:

    """
    perl /usr/local/bin/gsnapSam2Junctions.pl --is_bam --input_file ${bam} --output_file junctions.tab
    """
}

// work flow difination
workflow {
   // call the QC step
    fastqc = QualityControl(singleEndReads)

    index_ch = HisatIndex(params.reference)
   
   
   // Trimming based on paired or not
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(read_pairs_ch)

        hisat = HisatMappingPairedEnd(trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 

    }
    else 
    {
        trimm = SingleEndTrimming(singleEndReads,params.adapters)

        hisat = HisatMappingSingleEnd(trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
    }
    
    
    sortedByName = SortBams(hisat)

    if (params.isStranded) {

        hisatCount = HtseqCountingStranded(sortedByName)

    } else
    {
        hisatCount = HtseqCountingUnstranded(sortedByName)
    }

    spliceCounts = SpliceCrossingReads(hisat)
}
