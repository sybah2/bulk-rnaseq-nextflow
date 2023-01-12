#!/bin/bash
set -e
set -u

#java -jar trimmomatic.jar PE -trimlog ${workingDir}/trimLog $mateA $mateB -$mateAEncoding -baseout ${workingDir}/trimmedReads/${sampleName} /
#ILLUMINACLIP:\$GUS_HOME/data/DJob/DistribJobTasks/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20


# Single end trimming
# trimmomatic SE -phred33 ReadA_1.fastq.gz ReadA_1_trimmed.fastq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Paired end
#java -jar trimmomatic-0.39.jar  PE -trimlog Fastq/log Fastq/schizont_1.fq Fastq/schizont_2.fq Fastq/schizont_trim_paired_1.fq Fastq/schizont_trim_unpaired_1.fq Fastq/schizont_trim_paired_2.fq Fastq/schizont_trimm_unpaired_2.fq  ILLUMINACLIP:/Users/saikouybah/Documents/RNA_Seq_Nexflow/Adapters/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20


##java -jar $projectDir/trimmomatic-0.39.jar  PE -trimlog TrimLogs.txt ${reads[0]} ${reads[1]} ${sample_id}_trim_paired_1.fq ${sample_id}_trim_unpaired_1.fq ${sample_id}_trim_paired_2.fq ${sample_id}_trim_unpaired_2.fq  ILLUMINACLIP:$projectDir/Adapters/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20



### worked
#trimmomatic  PE -trimlog TrimLogs.txt ${reads[0]} ${reads[1]} ${sample_id}_trim_paired_1.fq ${sample_id}_trim_unpaired_1.fq ${sample_id}_trim_paired_2.fq ${sample_id}_trim_unpaired_2.fq  ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

# process QualityControl{

#     publishDir "${params.results}/QC_reports/", mode: 'copy'


#     input:
#     tuple val(sample_id), path(reads)

#     output:
#     path "${sample_id}_QC_report"

#     script:

#     """
#     FastQC.sh  "$reads", "${sample_id}"
#     """


# }


#java -jar $projectDir/trimmomatic-0.39.jar SE -trimlog ${sample_id}_Trimlogs.txt ${reads} ${sample_id}_trim.fq ILLUMINACLIP:${projectDir}/Adapters/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

hisat2 -x ${index} --max-intronlen ${intronLenght} -1 ${reads[1]} -2 ${reads[2]}  | samotools view - bS -  | samtools sort - > ${sample_id}_sorted.bam



