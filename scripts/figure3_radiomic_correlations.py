"""Figure 3: Spearman correlations among the 25 T2 radiomic features within NED and BCR patients.
Asterisks mark nominal p < 0.05 for the within-group correlation (not corrected; descriptive).
Usage: python figure3_radiomic_correlations.py <figures_dir>
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from paths import DATA
FIG = sys.argv[1]; BASE = DATA
feats = ['Maximum3DDiameter', 'Maximum2DDiameterSlice', 'MinorAxis', 'Volume', 'MajorAxis', 'LeastAxis', 'SurfaceArea',
         'Energy', 'RobustMeanAbsoluteDeviation', 'Maximum', 'RootMeanSquared', 'Minimum', 'StandardDeviation', 'Range',
         'Kurtosis', 'SumAverage', 'Autocorrelation', 'AverageIntensity', 'GrayLevelNonUniformity', 'LongRunEmphasis',
         'ShortRunHighGrayLevelEmphasis', 'RunLengthNonUniformity', 'LongRunHighGrayLevelEmphasis', 'RunEntropy',
         'HighGrayLevelRunEmphasis']
d = pd.read_csv(BASE + 'radio_T2v1.csv', decimal=',')
for c in d.columns[1:]: d[c] = pd.to_numeric(d[c].astype(str).str.replace(',', '.'), errors='coerce')
INK2 = '#52514e'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})
fig, axes = plt.subplots(1, 2, figsize=(17, 8.2), gridspec_kw={'wspace': 0.42})
for ax, (grp, name) in zip(axes, [(0, 'NED'), (1, 'BCR')]):
    g = d[d['Biochemical Recurrence'] == grp][feats]
    n = len(g)
    rho = g.corr(method='spearman').values
    p = np.array([[spearmanr(g[a], g[b], nan_policy='omit').pvalue if a != b else 1 for b in feats] for a in feats])
    im = ax.imshow(rho, cmap='RdBu_r', vmin=-1, vmax=1)
    for i in range(len(feats)):
        for j in range(len(feats)):
            if i != j and p[i, j] < 0.05:
                ax.text(j, i + 0.18, '*', ha='center', va='center', fontsize=8, color='#0b0b0b')
    ax.set_xticks(range(len(feats))); ax.set_yticks(range(len(feats)))
    ax.set_xticklabels(feats, rotation=90, fontsize=7.5); ax.set_yticklabels(feats, fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_title(f'{name} (n = {n})', fontsize=10, loc='left')
    ax.text(-0.02, 1.02, 'AB'[grp], transform=ax.transAxes, fontsize=14, fontweight='bold', ha='right')
cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02); cb.set_label('Spearman correlation', fontsize=9); cb.ax.tick_params(labelsize=8)
for ext in ('pdf', 'png'):
    fig.savefig(f'{FIG}/Figure3.{ext}', dpi=300, bbox_inches='tight')
print('Figure 3 written')
