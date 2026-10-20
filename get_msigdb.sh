#!/bin/bash
# Downloads the MSigDB v2023.2 gene sets used in the analysis (Hallmark and GO Biological Process, gene symbols).
set -euo pipefail
OUT=${1:-data/msigdb}; mkdir -p "$OUT"
URL=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs
for f in h.all.v2023.2.Hs.symbols.gmt c5.go.bp.v2023.2.Hs.symbols.gmt; do
  curl -fsSL "$URL/$f" -o "$OUT/$f"
done
echo "MSigDB gene sets saved in $OUT"
