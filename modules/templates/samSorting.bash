#!/usr/bin/env bash

set -euo pipefail

samtools sort ${sam} -o ${split_name}_sorted.bam