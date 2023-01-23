#!/bin/bash
set -e
set -u

#projectDir=${1}
encoding=${1}
read1=${2}
read2=${3}
read1_paired=${4}
read1_unpaired=${5}
read2_paired=${6}
read2_unpaired=${7}
adapters=${8}
sample_id=${9}


#java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE -trimlog ${sample_id} ${read1} ${read2} ${read1_paired} ${read2_unpaired} ${read2_paired} ${read2_unpaired} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20


#PairedEndTrimming.sh ${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz ${params.adaptersPE} ${sample_id}_Trimlog.txt 

java -jar /usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${encoding} -trimlog ${sample_id} ${read1} ${read2} ${read1_paired} ${read2_unpaired} ${read2_paired} ${read2_unpaired} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20