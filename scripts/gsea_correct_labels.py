"""Pre-ranked GSEA on the corrected DESeq2 results (BCR vs NED, verified labels).

Ranking metric: DESeq2 Wald statistic (BCR vs NED); duplicated symbols keep the largest |stat|.
Gene sets: MSigDB Hallmark 2023.2 (primary) and GO Biological Process 2023.2 (secondary), symbols.
Usage: python gsea_correct_labels.py <deseq2_csv> <gmt_dir> <out_dir>
"""
import sys, warnings
import pandas as pd, gseapy as gp
warnings.filterwarnings('ignore')
deseq, gmt, out = sys.argv[1:4]
d = pd.read_csv(deseq).dropna(subset=['stat', 'Gene name'])
d['abs'] = d['stat'].abs()
rnk = d.sort_values('abs', ascending=False).drop_duplicates('Gene name').set_index('Gene name')['stat'].sort_values(ascending=False)
for tag, f, mn, mx in [('hallmark', 'h.all.v2023.2.Hs.symbols.gmt', 15, 500), ('go_bp', 'c5.go.bp.v2023.2.Hs.symbols.gmt', 15, 500)]:
    res = gp.prerank(rnk=rnk, gene_sets=f'{gmt}/{f}', min_size=mn, max_size=mx, permutation_num=10000,
                     seed=2026, threads=16, outdir=None, verbose=False).res2d
    res = res[['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'Tag %', 'Lead_genes']].sort_values('NOM p-val')
    res.to_csv(f'{out}/gsea_{tag}.csv', index=False)
    q = pd.to_numeric(res['FDR q-val'])
    print(f'== {tag}: {len(res)} sets tested; FDR<0.25: {(q < 0.25).sum()}; FDR<0.05: {(q < 0.05).sum()}')
    print(res.head(12)[['Term', 'NES', 'NOM p-val', 'FDR q-val']].to_string(index=False))
