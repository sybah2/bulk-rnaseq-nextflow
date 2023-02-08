#!/usr/bin/env bash

set -euo pipefail

gsnapSam2Junctions.pl --is_bam --input_file ${bam} --output_file junctions.tab
