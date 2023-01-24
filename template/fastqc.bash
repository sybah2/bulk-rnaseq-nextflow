#!/usr/bin/env bash

set -euo pipefail

mkdir ${sample_id}_QC
fastqc ${params.fastqcExtension}  -o "${sample_id}_QC" --extract 