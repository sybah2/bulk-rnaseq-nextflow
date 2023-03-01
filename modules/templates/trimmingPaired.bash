#!/usr/bin/env bash

set -euo pipefail

mateAEncoding=\$(cat ${quality_check_out})
java -jar $projectDir/trimmomatic-0.39.jar PE -\$mateAEncoding -trimlog ${sample_id}_Trimlog.txt \
${reads[0]} ${reads[1]} ${sample_id}_paired_1.fq.gz  ${sample_id}_unpaired_1.fq.gz ${sample_id}_paired_2.fq.gz ${sample_id}_unpaired_2.fq.gz \
ILLUMINACLIP:$projectDir/data/Adapters/All_adaptors-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

#/usr/local/Trimmomatic-0.39
#/usr/local/bin/