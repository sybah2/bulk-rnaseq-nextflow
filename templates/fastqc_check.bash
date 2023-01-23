#!/usr/bin/env bash

set -euo pipefail

perl /usr/local/bin/fastqc_check.pl ${fastqc_out} quality_check_out