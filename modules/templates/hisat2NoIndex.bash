#!/usr/bin/env bash

set -euo pipefail
TMP=$params.hisat2Index
FILES=\$TMP*
for f in \$FILES; do cp "\$f" "\${f#\$TMP}" ; done
