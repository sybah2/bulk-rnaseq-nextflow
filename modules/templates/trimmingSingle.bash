#!/usr/bin/env bash

set -euo pipefail

mateAEncoding=\$(cat ${quality_check_out})
java -jar/usr/local/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads ${task.cpus} -\$mateAEncoding -trimlog ${sample_id}_Trimlog.txt \
${reads} ${sample_id}_trim.fq.gz ILLUMINACLIP:/usr/local/bin/All_adaptors-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

