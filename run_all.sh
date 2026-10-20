#!/bin/bash
# Reproduces the results, tables and figures of the article from the input tables in data/ (see data/README.md).
#   bash run_all.sh                              # reads ./data, writes ./results and ./figures
#   RGBCR_DATA=/path/to/tables bash run_all.sh   # input tables in another folder
#   PYTHON=/path/to/python bash run_all.sh       # a specific interpreter
set -euo pipefail
ROOT=$(cd "$(dirname "$0")" && pwd)
PY=${PYTHON:-python3}
C=$ROOT/scripts; R=$ROOT/results; F=$ROOT/figures; G=$ROOT/data/msigdb
mkdir -p "$R/explanatory" "$F"
[ -f "$G/h.all.v2023.2.Hs.symbols.gmt" ] || bash "$ROOT/get_msigdb.sh" "$G"
cd "$R"
$PY "$C/deseq_correct.py" "$R/deseq2_BCR_vs_NED_correct_labels.csv"                 # differential expression (Figure 2A, Table S1)
$PY "$C/gsea_correct_labels.py" "$R/deseq2_BCR_vs_NED_correct_labels.csv" "$G" "$R/explanatory"   # pre-ranked GSEA (Table S3B)
$PY "$C/gsea_phenotype_perm.py" "$G" "$R/explanatory"                                 # GSEA, phenotype permutation (Figure 2B, Table S3A)
$PY "$C/radiomics_fdr.py" "$R"                                                        # radiomics with FDR (Table 2, Table S2)
$PY "$C/performance_table.py" "$R"                                                    # cross-validated performance (Table S5)
$PY "$C/deseq_topk_selection.py" "$R"                                                 # gene pre-selection bias (response to reviewers)
$PY "$C/top20_explanatory.py" "$R"                                                    # re-selection frequency of the top-20 genes (Table S4)
$PY "$C/figure1_workflow.py" "$F"
$PY "$C/figure2_de_gsea.py" "$R" "$F"
$PY "$C/figure3_radiomic_correlations.py" "$F"
$PY "$C/figure4_explanatory.py" "$R" "$F"
echo "done: results in $R, figures in $F"
