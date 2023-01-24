#!/usr/bin/env bash

set -euo pipefail

samtools sort -n ${bam} -o "${sample_id}_sortedByName.bam"