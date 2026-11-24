#!/usr/bin/env bash
# Run the whole pipeline end to end. Usage: ./run_pipeline.sh [first_step] [last_step]
set -euo pipefail
cd "$(dirname "$0")"
first=${1:-0}; last=${2:-6}
for script in scripts/0[0-6]_*.py; do
  n=$(basename "$script" | cut -c1-2); n=$((10#$n))
  if (( n >= first && n <= last )); then
    python "$script"
  fi
done
