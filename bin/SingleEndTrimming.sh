#!/bin/bash
set -e
set -u

projectDir=${1}
read=${2}
read_trim=${3}
adapters=${4}
sample_id=${5}


java -jar $projectDir/trimmomatic-0.39.jar SE -trimlog ${sample_id} ${read} ${read_trim} ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20