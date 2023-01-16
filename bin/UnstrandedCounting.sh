#!/bin/bash
set -e
set -u

bam=${1}
annotation=${2}
output1=${3}
output2=${4}

htseq-count -a 0 -f bam -s no -t exon -i gene_id ${bam} ${annotation} > ${output1}

htseq-count -a 0 -f bam --nonunique all -s no -t exon -i gene_id ${bam} ${annotation} > ${output2}
