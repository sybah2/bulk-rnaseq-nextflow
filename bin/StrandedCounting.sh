#!/bin/bash
set -e
set -u

bam=${1}
annotation=${2}
output1=${3}
output2=${4}
output3=${5}
output4=${6}

htseq-count -a 0 -f bam -s reverse -t exon -i gene_id ${bam} ${annotation} > ${output1}

htseq-count -a 0 -f bam -s yes -t exon -i gene_id ${bam} ${annotation} > ${output2}

htseq-count -a 0 -f bam --nonunique all -s reverse -t exon -i gene_id ${bam} ${annotation} > ${output3}

htseq-count -a 0 -f bam --nonunique all -s yes -t exon -i gene_id ${bam} ${annotation} > ${output4}