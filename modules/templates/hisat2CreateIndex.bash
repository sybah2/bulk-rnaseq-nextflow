#!/usr/bin/env bash

set -euo pipefail

if [ "$createIndex" = true ]; then

    hisat2-build ${reference} genomeIndex
    
else

    TMP=$hisat2Index
    FILES=\$TMP/*
    for f in \$FILES; do cp "\$f" "." ; done

fi
