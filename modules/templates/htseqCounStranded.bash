#!/usr/bin/env bash

set -euo pipefail

if ["isCds" == true]; then

    htseq-count -a 0 -f bam -s reverse -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.firststrand.counts"

    htseq-count -a 0 -f bam -s yes -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.secondstrand.counts"

    htseq-count -a 0 -f bam --nonunique all -s reverse -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.firststrand.nonunique.counts" 

    htseq-count -a 0 -f bam --nonunique all -s yes -t CDS -i gene_id ${bam} ${annotation} > "genes.htseq-union.secondstrand.nonunique.counts"

else 
    htseq-count -a 0 -f bam -s reverse -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.firststrand.counts"

    htseq-count -a 0 -f bam -s yes -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.secondstrand.counts"

    htseq-count -a 0 -f bam --nonunique all -s reverse -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.firststrand.nonunique.counts" 

    htseq-count -a 0 -f bam --nonunique all -s yes -t exon -i gene_id ${bam} ${annotation} > "genes.htseq-union.secondstrand.nonunique.counts"

fi