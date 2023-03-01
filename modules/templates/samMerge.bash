#!/usr/bin/env bash

set -euo pipefail

samtools merge -f  ${sampleID}_temp.bam *.bam

samtools view -bS ${sampleID}_temp.bam | samtools sort - > ${sampleID}.bam

rm ${sampleID}_temp.bam

