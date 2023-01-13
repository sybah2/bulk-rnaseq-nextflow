#!/bin/bash
set -e
set -u

index=${1}
intronLenght=${2}
read1=${3}
read2=${4}
sample_id=${5}


hisat2 -x ${index} --max-intronlen ${intronLenght} -1 ${read1} -2 ${read2} 2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}_sorted.bam