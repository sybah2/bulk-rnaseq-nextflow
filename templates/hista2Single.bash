#!/usr/bin/env bash

set -euo pipefail

hisat2 -x ${index} --max-intronlen ${params.intronLenght} -U ${paired1}  2>hisat2.log | samtools view -bS - | samtools sort - > ${sample_id}_sorted.bam