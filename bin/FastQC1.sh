#!/bin/bash
set -e
set -u

reads=${1}
sample_id=${2}


mkdir ${sample_id}_QC_report
fastqc ${reads} -o ${sample_id}_QC_report --extract
