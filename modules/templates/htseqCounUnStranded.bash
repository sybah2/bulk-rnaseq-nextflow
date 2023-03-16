#!/usr/bin/env bash

set -euo pipefail

if ["isCds" == true]; then
    htseq-count -a 0 -f bam -s no -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.unstranded.counts" 

    htseq-count -a 0 -f bam --nonunique all -s no -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.unstranded.nonunique.counts"

else

    htseq-count -a 0 -f bam -s no -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.unstranded.counts" 

    htseq-count -a 0 -f bam --nonunique all -s no -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.unstranded.nonunique.counts"

fi