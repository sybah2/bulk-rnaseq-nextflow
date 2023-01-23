#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the channel for single and paired end reads
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

    """
    hisat2-build ${reference} genomeIndex
    """  
}

// FastQC quality control process defination
process QualityControl {

    tag {sample_id}

    container = 'biocontainers/fastqc:v0.11.9_cv7'

    //publishDir "${params.results}/${sample_id}", mode: 'copy'
    input:
    path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_QC")

    script:
    sample_id = reads.getSimpleName()

    """
    mkdir ${sample_id}_QC
    fastqc -o ${sample_id}_QC -f fastq -q ${reads} --extract 
    
    """
}

process Fastqc_check {

    tag {sample_id}

    container = 'veupathdb/shortreadaligner'
    

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastqc_out)

    output:
    path("quality_check_out")

    script:

    """
    perl /usr/local/bin/fastqc_check.pl ${fastqc_out} "quality_check_out"
    """

}

// Paired end trimming (with Trimmomatic) process defination 
process PaireEndTrimming {

    //container = 'veupathdb/shortreadaligner'

    tag {sample_id}

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    path("quality_check_out")
    tuple val(sample_id), path(reads)
    
    

    output:
    tuple val(sample_id), path("${sample_id}_paired_1.fq.gz"), path("${sample_id}_paired_2.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log
    
    script:
    """  
    mateAEncoding=\$(cat "quality_check_out")
    java -jar ${projectDir}/trimmomatic-0.39.jar PE -phred64 -trimlog ${sample_id}_Trimlog.txt ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz ILLUMINACLIP:${params.adaptersPE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    """
}

// Single end trimming (With Trimmomatic) process defination
process SingleEndTrimming {

    //container = 'veupathdb/shortreadaligner'

    tag {sample_id}
    
    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    path("quality_check_out")
    path(reads)
    

    output:
    tuple val(sample_id), path("${sample_id}_trim.fq.gz"), emit: trimmed_fastqs
    path("${sample_id}_Trimlog.txt"), emit: trimm_log

    script:
    sample_id = reads.getSimpleName()

    """
    mateAEncoding=\$(cat "quality_check_out")
    java -jar ${projectDir}/trimmomatic-0.39.jar SE -\$mateAEncoding -trimlog ${sample_id}_Trimlog.txt ${reads} ${sample_id}_trim.fq.gz ILLUMINACLIP:${params.adaptersSE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
    
    """

}

process HisatMappingPairedEnd{
    container = 'veupathdb/shortreadaligner'

    //publishDir "${params.results}/${sample_id}", mode: 'copy'
    tag {sample_id}

    input:
    tuple val(sample_id), path(paired1), path(paired2)
    val index
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    
    """
    hisat2 -x ${index} --max-intronlen ${params.intronLenght} -1 ${paired1} -2 ${paired2} 2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}.bam
    """

}

process HisatMappingSingleEnd{
    tag {sample_id}
    container = 'veupathdb/shortreadaligner'

    //publishDir "${params.results}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(paired1)
    val index 
    path 'genomeIndex.*.ht2'

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    hisat2 -x ${index} --max-intronlen ${params.intronLenght} -U ${paired1}  2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}.bam
    
    """

}

process SortBams{
    tag {sample_id}
    container = 'veupathdb/shortreadaligner'
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
    
    container = 'genomicpariscentre/htseq'
    tag {sample_id}
    
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
    
    container = 'genomicpariscentre/htseq'
    tag {sample_id}
    
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
    
    tag {sample_id}

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

process BedFileStatsStranded {
    tag {sample_id}
    
    publishDir "${params.results}/${sample_id}", mode: 'copy'

    container = 'saikou'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple path("*.bed"), path("mappingStats.txt")
    
    script:

    """
    perl /usr/local/bin/gsnapSplitBam.pl --mainResultDir . --strandSpecific 1 --isPairedEnd 0 --bamFile ${bam}

    """
}

// work flow difination
workflow {
   // call the QC step
    fastqc = QualityControl(singleEndReads)
    //fastqc.view()
    chech_fastq = Fastqc_check(fastqc)
    //chech_fastq.view()
    
    index_ch = HisatIndex(params.reference)
   
   
   // Trimming based on paired or not
    if (params.isPaired) 
    {
        trimm = PaireEndTrimming(chech_fastq,read_pairs_ch)

        hisat = HisatMappingPairedEnd(trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
        hisat.view()

    }
    else 
    {
        trimm = SingleEndTrimming(chech_fastq,singleEndReads)

        hisat = HisatMappingSingleEnd(trimm.trimmed_fastqs, index_ch.genome_index_name, index_ch.ht2_files) 
        hisat.view()
    }
    
  
    sortedByName = SortBams(hisat)

    if (params.isStranded) {

        hisatCount = HtseqCountingStranded(sortedByName)

    } else
    {
        hisatCount = HtseqCountingUnstranded(sortedByName)
    }

    spliceCounts = SpliceCrossingReads(hisat)

    befFileStats = BedFileStatsStranded(hisat)


}
