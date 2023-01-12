#!/bin/bash
set -e
set -u

projectDir=${1}
read1=${2}
read2=${3}
read1_paired=${4}
read2_paired=${5}
adapters=${6}
sample_id=${7}


java -jar ${projectDir}/trimmomatic-0.39.jar  PE -trimlog ${sample_id} ${read1} ${read2} ${read1_paired} ${read2_paired}  ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20