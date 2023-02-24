#!/usr/bin/env bash

set -euo pipefail

fasterq-dump --gzip --split-3 ${sra}
