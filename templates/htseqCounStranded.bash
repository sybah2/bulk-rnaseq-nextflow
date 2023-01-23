#!/usr/bin/env bash

set -euo pipefail

##StrandedCounting.sh ${bam} ${params.annotation} "genes.htseq-union.firststrand.counts" "genes.htseq-union.secondstrand.counts" "genes.htseq-union.firststrand.nonunique.counts" "genes.htseq-union.secondstrand.nonunique.counts"


htseq-count -a 0 -f bam -s reverse -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.firststrand.counts"

htseq-count -a 0 -f bam -s yes -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.secondstrand.counts"

htseq-count -a 0 -f bam --nonunique all -s reverse -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.firststrand.nonunique.counts" 

htseq-count -a 0 -f bam --nonunique all -s yes -t exon -i gene_id ${bam} ${params.annotation} > "genes.htseq-union.secondstrand.nonunique.counts"