"""GSEA with phenotype (label) permutation: more conservative than pre-ranked gene-set permutation
when genes are correlated. Hallmark 2023.2; metric signal-to-noise; log2 DESeq2-normalised expression.
Usage: python gsea_phenotype_perm.py <gmt_dir> <out_dir>"""
import sys, warnings, numpy as np, pandas as pd, gseapy as gp
warnings.filterwarnings('ignore')
from paths import DATA
gmt, out = sys.argv[1:3]; BASE = DATA
cnt = pd.read_csv(BASE + 'genes_with_ensembl_ids.csv', low_memory=False)
cnt = cnt.drop(columns=[c for c in cnt.columns if c in ('Ensembl ID', 'Unnamed: 36')]).drop_duplicates('Gene name').set_index('Gene name').apply(pd.to_numeric, errors='coerce').fillna(0)  # one copy per gene: BioMart merge duplicated rows
cnt.columns = [c.upper() for c in cnt.columns]; cnt = cnt[((cnt >= 10).mean(axis=1) >= 0.10)]
lg = np.log(cnt.replace(0, np.nan)); ok = lg.notna().all(axis=1)
sf = np.exp(lg[ok].sub(lg[ok].mean(axis=1), axis=0).median()); E = np.log2(cnt / sf + 1)
lab = pd.read_csv(BASE + 'radio_T2.csv', decimal=','); lab.index = lab['Samples'].str.replace('-', '').str.upper()
cls = ['BCR' if lab.loc[s, 'Biochemical Recurrence'] == 1 else 'NED' for s in E.columns]
res = gp.gsea(data=E, gene_sets=f'{gmt}/h.all.v2023.2.Hs.symbols.gmt', cls=cls, permutation_type='phenotype',
              permutation_num=5000, method='signal_to_noise', min_size=15, max_size=500, seed=2026, threads=16, outdir=None).res2d
res = res[['Term', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes']].sort_values('NOM p-val')
res.to_csv(f'{out}/gsea_hallmark_phenotype_perm.csv', index=False)
q = pd.to_numeric(res['FDR q-val'])
print(f'phenotype permutation: FDR<0.25: {(q<0.25).sum()}  FDR<0.05: {(q<0.05).sum()}')
print(res.head(12)[['Term', 'NES', 'NOM p-val', 'FDR q-val']].to_string(index=False))
