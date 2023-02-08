#!/usr/bin/env bash

set -euo pipefail

mkdir "${sample_id}_fastqc_out"
fastqc -o "${sample_id}_fastqc_out" ${reads} --extract