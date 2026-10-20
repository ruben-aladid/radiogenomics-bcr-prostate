# Input data

The processed input tables are not distributed with the code. Raw sequencing reads for all 42 patients are available in the NCBI Sequence Read Archive under BioProject [PRJNA1479297](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1479297). For access to the processed tables, see the Data Availability Statement of the article.

Place the three tables below in this folder, or point `RGBCR_DATA` to the folder that contains them. The MSigDB gene sets are downloaded to `data/msigdb/` by `get_msigdb.sh` the first time `run_all.sh` runs.

| File | Content | Format |
|---|---|---|
| `genes_with_ensembl_ids.csv` | Gene-level read counts of the index lesion (Ion AmpliSeq Transcriptome, one amplicon per gene), annotated with Ensembl BioMart | Comma-separated. Columns: `Ensembl ID`, `Gene name`, then one column of integer counts per patient. Genes duplicated by the annotation are collapsed to one row per gene symbol by the scripts. |
| `radio_T2.csv` | Clinical outcome per patient | Comma-separated, decimal comma. Columns used: `Samples` (patient identifier; hyphens are removed and the identifier is upper-cased to match the count table) and `Biochemical Recurrence` (1 = BCR, 0 = NED). Other columns are ignored. |
| `radio_T2v1.csv` | T2-weighted PyRadiomics features of the index lesion (94 shape, first-order and texture features) | Comma-separated, decimal comma. Columns: `Samples`, `Biochemical Recurrence`, then one column per feature. |

Patient identifiers must match across the three tables (the scripts map `P1755` in the radiomic table to `P17551` in the count table).
