#!/usr/bin/env bash

set -euo pipefail

TMP=${params.hisat2Index}
FILES=\$TMP/*ht2
for f in \$FILES; do cp "\$f" "." ; done