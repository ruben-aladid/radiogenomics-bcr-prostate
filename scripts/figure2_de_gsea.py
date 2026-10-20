"""Figure 2 (differential expression + Hallmark GSEA) generated
directly from the analysis outputs. Writes PDF (vector) and PNG (300 dpi).
Usage: python figure2_de_gsea.py <results_dir> <figures_dir>
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shap

RC = sys.argv[1]; EX = f'{RC}/explanatory'; FIG = sys.argv[2]
UP, DOWN, NS = '#e34948', '#2a78d6', '#b5b4ae'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e0'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False,
                     'pdf.fonttype': 42, 'svg.fonttype': 'none'})

def panel(ax, letter):
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=13, fontweight='bold', color=INK, ha='right', va='bottom')

# ---------------- Figure 2
de = pd.read_csv(f'{RC}/deseq2_BCR_vs_NED_correct_labels.csv').dropna(subset=['padj', 'log2FoldChange'])
de['nlp'] = -np.log10(de['padj'])
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.9})
sig = de['padj'] < 0.1
col = np.where(~sig, NS, np.where(de['log2FoldChange'] > 0, UP, DOWN))
a.scatter(de['log2FoldChange'], de['nlp'], c=col, s=np.where(sig, 14, 5), linewidths=0, alpha=0.9, rasterized=True)
xmin, xmax = de['log2FoldChange'].min(), de['log2FoldChange'].max()
a.set_xlim(xmin - 2.6, xmax + 2.6)
for y, lbl in [(-np.log10(0.05), 'FDR 0.05'), (-np.log10(0.1), 'FDR 0.10')]:
    a.axhline(y, color=INK2, lw=0.8, ls='--')
    a.text(0, y + 0.015, lbl, fontsize=7.5, color=INK2, va='bottom', ha='center')
# label the 6 genes with the largest |log2FC| among FDR < 0.10 on each side, stacked in the margins
ymax = de['nlp'].max()
for side, sel in [(1, de[sig & (de['log2FoldChange'] > 0)]), (-1, de[sig & (de['log2FoldChange'] < 0)])]:
    sel = sel.reindex(sel['log2FoldChange'].abs().sort_values(ascending=False).index).head(6).sort_values('nlp', ascending=False)
    for i, (_, r) in enumerate(sel.iterrows()):
        tx = xmax + 0.6 if side > 0 else xmin - 0.6
        a.annotate(r['Gene name'], (r['log2FoldChange'], r['nlp']), xytext=(tx, ymax - 0.11 * i),
                   fontsize=7.5, color=INK, ha='left' if side > 0 else 'right', va='center',
                   arrowprops=dict(arrowstyle='-', color=INK2, lw=0.5, shrinkA=2, shrinkB=2))
a.set_xlabel('log2 fold change (BCR vs NED)'); a.set_ylabel('-log10 adjusted p-value')
a.grid(color=GRID, lw=0.6); a.set_axisbelow(True)
n_up, n_dn = int((sig & (de['log2FoldChange'] > 0)).sum()), int((sig & (de['log2FoldChange'] < 0)).sum())
a.legend(handles=[plt.Line2D([], [], marker='o', ls='', color=UP, label=f'Higher in BCR, FDR < 0.10 (n = {n_up})'),
                  plt.Line2D([], [], marker='o', ls='', color=DOWN, label=f'Lower in BCR, FDR < 0.10 (n = {n_dn})'),
                  plt.Line2D([], [], marker='o', ls='', color=NS, label='Not significant')],
         fontsize=7.5, frameon=False, loc='upper left', bbox_to_anchor=(0, -0.13), ncol=1)
panel(a, 'A')

g = pd.read_csv(f'{EX}/gsea_hallmark_phenotype_perm.csv')
g['NOM p-val'] = pd.to_numeric(g['NOM p-val']); g['FDR q-val'] = pd.to_numeric(g['FDR q-val'])
g = g[g['NOM p-val'] < 0.10].sort_values('NES')
ACR = {'Uv': 'UV', 'Kras': 'KRAS', 'Tnfa': 'TNFα', 'Nfkb': 'NF-κB', 'Il2': 'IL2', 'Stat5': 'STAT5', 'Myc': 'MYC',
       'E2f': 'E2F', 'G2m': 'G2M', 'Mtorc1': 'mTORC1', 'Pi3k': 'PI3K', 'Akt': 'AKT', 'Mtor': 'mTOR', 'Tgf': 'TGF',
       'Dna': 'DNA', 'Il6': 'IL6', 'Jak': 'JAK', 'Stat3': 'STAT3', 'P53': 'p53', 'Wnt': 'WNT', 'Notch': 'NOTCH'}
g['label'] = g['Term'].str.replace('HALLMARK_', '').str.replace('_', ' ').str.lower().str.capitalize().apply(
    lambda s: ' '.join(ACR.get(w.capitalize(), w) for w in s.split(' ')))
b.barh(g['label'], g['NES'], color=np.where(g['NES'] > 0, UP, DOWN), height=0.62)
for yi, (_, r) in enumerate(g.iterrows()):
    b.text(r['NES'] - 0.05 if r['NES'] < 0 else r['NES'] + 0.05, yi, f"q = {r['FDR q-val']:.2f}", va='center',
           ha='right' if r['NES'] < 0 else 'left', fontsize=7.5, color=INK2)
b.axvline(0, color=INK2, lw=0.8)
b.set_xlim(min(-2.4, g['NES'].min() - 0.7), max(1.0, g['NES'].max() + 0.7))
b.set_xlabel('Normalised enrichment score (phenotype permutation)')
b.tick_params(axis='y', length=0, labelsize=8); b.spines['left'].set_visible(False)
b.grid(axis='x', color=GRID, lw=0.6); b.set_axisbelow(True)
b.set_title('Hallmark gene sets, nominal p < 0.10', fontsize=9, color=INK, loc='left')
panel(b, 'B')
for ext in ('pdf', 'png'):
    fig.savefig(f'{FIG}/Figure2.{ext}', dpi=300, bbox_inches='tight')
plt.close(fig)

print('Figure 2 written')
