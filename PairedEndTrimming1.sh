#!/bin/bash
set -e
set -u

projectDir=${1}
read1=${2}
read2=${3}
read1_paired=${4}
read1_unpaired=${5}
read2_paired=${6}
read2_unpaired=${7}
adapters=${8}
sample_id=${9}


java -jar ${projectDir}/trimmomatic-0.39.jar  PE -trimlog ${sample_id} ${read1} ${read2} ${read1_paired} ${read2_unpaired} ${read2_paired} ${read2_unpaired} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
