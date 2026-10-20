# Radiogenomic analysis of biochemical recurrence in localized prostate cancer

Analysis code for the article *Machine Learning-Based Analysis of Biochemical Recurrence in Localized Prostate Cancer Using Radiogenomic Data: An Exploratory Study* (Ajami T, et al.; under review in *Biomolecules*). The citation will be added on publication.

The study profiles the transcriptome (Ion AmpliSeq Transcriptome) and the T2-weighted MRI radiomic features of the same MRI-defined index lesion in patients treated with radical prostatectomy, and explores their association with biochemical recurrence (BCR). The analysed cohort comprises 35 patients with paired data (11 BCR, 24 NED).

## What the code does

- **Differential expression** between BCR and NED with PyDESeq2, and gene set enrichment with GSEApy (MSigDB Hallmark; phenotype permutation as the primary analysis, pre-ranked as secondary).
- **Radiomics:** between-group comparisons (Mann–Whitney U) and differences in Spearman correlations (Fisher's z), with Benjamini–Hochberg correction per family of tests.
- **Explanatory models:** XGBoost fitted on all patients with SHAP values. These models are descriptive and are not used to estimate performance.
- **Predictive performance:** 5-fold stratified cross-validation repeated 10 times, with gene selection (DESeq2), imputation and scaling performed inside the training folds only.
- **Selection bias:** the same models with genes selected on the full cohort before cross-validation, with true and permuted outcomes, to show how much of a high AUC comes from information leakage.

## Scripts and outputs

| Script | Output (in `results/` or `figures/`) | Article |
|---|---|---|
| `deseq_correct.py` | `deseq2_BCR_vs_NED_correct_labels.csv` | Figure 2A, Table S1, Data S1 |
| `gsea_phenotype_perm.py` | `explanatory/gsea_hallmark_phenotype_perm.csv` | Figure 2B, Table S3A |
| `gsea_correct_labels.py` | `explanatory/gsea_hallmark.csv`, `explanatory/gsea_go_bp.csv` | Table S3B, Data S3–S4 |
| `radiomics_fdr.py` | `radiomics_mannwhitney_fdr.csv`, `radiomics_fisher_z_fdr.csv` | Table 2, Table S2, Data S2 |
| `performance_table.py` | `table3_cv_performance.csv`, `table3_gene_selection_frequency.csv` | Table S5, Data S5 |
| `top20_explanatory.py` | `explanatory/top20_explanatory_table.csv` | Table S4 |
| `deseq_topk_selection.py` | `deseq_topk_auc_summary.csv`, `deseq_topk_auc_all.csv`, `deseq_topk_gene_stability.csv` | Response to reviewers (selection bias) |
| `figure1_workflow.py` … `figure4_explanatory.py` | `Figure1`–`Figure4` (PDF and PNG); `figure4_*_shap.csv` | Figures 1–4 |

## Requirements

Python 3.11 and the packages in `requirements.txt`:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Data

Raw sequencing reads are available in the NCBI Sequence Read Archive under BioProject [PRJNA1479297](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479297). The processed input tables (gene counts, outcome and radiomic features) are not included in this repository; [`data/README.md`](data/README.md) describes the files and their format.

## Running

```bash
bash run_all.sh                              # input tables in ./data
RGBCR_DATA=/path/to/tables bash run_all.sh   # input tables elsewhere
```

`run_all.sh` downloads the MSigDB v2023.2 gene sets on first use (`get_msigdb.sh`), runs every analysis in order and writes the results to `results/` and the figures to `figures/`. A full run takes about 7 minutes on a 48-core workstation; most of the time is spent re-running DESeq2 inside the cross-validation folds.

Every stochastic step (cross-validation splits, bootstrap refits, permutations, XGBoost) uses a fixed seed. Small numerical differences may appear with other package versions.

## License

MIT (see `LICENSE`).
