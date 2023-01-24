#!/usr/bin/env bash

set -euo pipefail

#UnstrandedCounting.sh ${bam} ${params.annotation} "genes.htseq-union.unstranded.counts" "genes.htseq-union.unstranded.nonunique.counts"


htseq-count -a 0 -f bam -s no -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.unstranded.counts" 

htseq-count -a 0 -f bam --nonunique all -s no -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.unstranded.nonunique.counts"