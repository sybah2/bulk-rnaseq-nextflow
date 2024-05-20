#!/usr/bin/env bash

set -euo pipefail

gsnapSplitBam.pl --mainResultDir . \
                --strandSpecific ${strandSpecific} \
                --isPairedEnd ${isPairedEnd} \
                --bamFile ${bam}