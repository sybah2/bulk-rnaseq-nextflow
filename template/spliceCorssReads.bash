#!/usr/bin/env bash

set -euo pipefail

perl /usr/local/bin/gsnapSam2Junctions.pl --is_bam --input_file ${bam} --output_file junctions.tab