#!/usr/bin/env bash

set -euo pipefail

if [ "$createIndex" = true ]; then

    hisat2-build ${reference} ${organismAbbv}
    
else

    TMP=$hisat2Index
    FILES=\$TMP/*ht2
    for f in \$FILES; do cp "\$f" "." ; done

fi
