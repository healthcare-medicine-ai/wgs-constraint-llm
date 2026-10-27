#!/usr/bin/env bash
# Run every stage in order. Each stage reads the outputs of the ones before it.
# The whole run takes about four hours on one node with 128 GB of memory.
set -euo pipefail
cd "$(dirname "$0")"

for stage in pipelines/0*.py; do
    echo "== $stage"
    python "$stage"
done
