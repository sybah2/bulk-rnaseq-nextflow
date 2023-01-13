#!/bin/bash
set -e
set -u

index=${1}
intronLenght=${2}
read1=${3}
sample_id=${4}


hisat2 -x ${index} --max-intronlen ${intronLenght} -U ${read1}  2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}_sorted.bam