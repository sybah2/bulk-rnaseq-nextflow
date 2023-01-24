#!/usr/bin/env bash

set -euo pipefail

perl /usr/local/bin/gsnapSplitBam.pl --mainResultDir . --strandSpecific 0 --isPairedEnd 0 --bamFile ${bam}