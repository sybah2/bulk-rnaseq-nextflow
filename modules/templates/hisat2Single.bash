#!/usr/bin/env bash

set -euo pipefail
mateAEncoding=\$(cat ${quality_check_out})
hisat2 --\$mateAEncoding -x ${index} --max-intronlen ${params.intronLength} -U ${read}  2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}.bam