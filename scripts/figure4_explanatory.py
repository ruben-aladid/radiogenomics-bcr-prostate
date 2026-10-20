"""Figure 4: descriptive (non-predictive) explanatory models with SHAP, verified BCR labels.

All three models use the same XGBoost settings (original RNA model) and are fitted on all patients:
  A/B transcriptomic: 20 genes with smallest full-cohort DESeq2 adjusted p-values (n = 35)
  C/D radiomic:       94 T2 radiomic features, mean-imputed (n = 34)
  E/F combined:       the 20 genes + 94 radiomic features (n = 34)
SHAP values are computed on all patients. Stability: % of 200 bootstrap refits in which a feature is
among the model's 10 features with the largest mean |SHAP| (features with zero SHAP are not counted).
No performance metric is derived from these models (see Table 3 for cross-validated performance).
Usage: python figure4_explanatory.py <results_dir> <figures_dir>
"""
import sys, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shap
from xgboost import XGBClassifier
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')
from paths import DATA
RES, FIG = sys.argv[1], sys.argv[2]; BASE = DATA
rng = np.random.default_rng(2026)

cnt = pd.read_csv(BASE + 'genes_with_ensembl_ids.csv', low_memory=False)
cnt = cnt.drop(columns=[c for c in cnt.columns if c in ('Ensembl ID', 'Unnamed: 36')]).drop_duplicates('Gene name') \
         .set_index('Gene name').apply(pd.to_numeric, errors='coerce').fillna(0)
cnt.columns = [c.upper() for c in cnt.columns]; cnt = cnt[((cnt >= 10).mean(axis=1) >= 0.10)]
lg = np.log(cnt.replace(0, np.nan)); ok = lg.notna().all(axis=1)
E = np.log2(cnt / np.exp(lg[ok].sub(lg[ok].mean(axis=1), axis=0).median()) + 1).T
de = pd.read_csv(f'{RES}/deseq2_BCR_vs_NED_correct_labels.csv').dropna(subset=['padj'])
de['abs'] = de['stat'].abs()
top20 = de.sort_values(['padj', 'abs'], ascending=[True, False])['Gene name'].head(20).tolist()
r = pd.read_csv(BASE + 'radio_T2v1.csv', decimal=',')
r.index = r['Samples'].str.upper().replace({'P1755': 'P17551'})
R = r.drop(columns=['Samples', 'Biochemical Recurrence']).apply(
    lambda s: pd.to_numeric(s.astype(str).str.replace(',', '.'), errors='coerce')).dropna(axis=1, how='all')
R = pd.DataFrame(SimpleImputer(strategy='mean').fit_transform(R), index=R.index, columns=R.columns)
lab = pd.read_csv(BASE + 'radio_T2.csv', decimal=','); lab.index = lab['Samples'].str.replace('-', '').str.upper()
truth = lab['Biochemical Recurrence'].astype(int)
std = lambda X: pd.DataFrame(StandardScaler().fit_transform(X), index=X.index, columns=X.columns)

ids_c = [i for i in R.index if i in E.index]
blocks = {'transcriptomic': std(E[top20]), 'radiomic': std(R), 'combined': std(pd.concat([E.loc[ids_c, top20], R.loc[ids_c]], axis=1))}
P = dict(n_estimators=200, max_depth=5, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8, eval_metric='logloss', random_state=42)

out = {}
for key, X in blocks.items():
    y = truth.reindex(X.index).astype(int)
    sv = shap.TreeExplainer(XGBClassifier(**P).fit(X, y)).shap_values(X)
    imp = pd.Series(np.abs(sv).mean(axis=0), index=X.columns)
    freq = pd.Series(0.0, index=X.columns)
    for _ in range(200):
        idx = rng.choice(len(X), len(X), replace=True)
        if y.iloc[idx].nunique() < 2: continue
        sb = np.abs(shap.TreeExplainer(XGBClassifier(**P).fit(X.iloc[idx], y.iloc[idx])).shap_values(X.iloc[idx])).mean(axis=0)
        t = np.argsort(sb)[::-1][:10]; freq[X.columns[t[sb[t] > 0]]] += 1
    freq /= 200
    tab = pd.DataFrame({'mean_abs_shap': imp, 'pct_bootstrap_top10': (freq * 100).round(0)}).sort_values('mean_abs_shap', ascending=False)
    tab.to_csv(f'{RES}/figure4_{key}_shap.csv')
    out[key] = (sv, X, tab, len(y), int(y.sum()))
    print(key, 'n', len(y), 'BCR', int(y.sum()), '| top 8:', ', '.join(f'{g} ({p:.0f}%)' for g, p in tab['pct_bootstrap_top10'].head(8).items()))

INK, INK2, GRID, BLUE = '#0b0b0b', '#52514e', '#e6e5e0', '#2a78d6'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'pdf.fonttype': 42})
titles = {'transcriptomic': 'Transcriptomic: top-20 DESeq2 genes', 'radiomic': 'Radiomic: T2 features', 'combined': 'Combined: top-20 genes + T2 features'}
fig = plt.figure(figsize=(12.5, 14)); gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.25], hspace=0.42, wspace=0.75)
for row, key in enumerate(['transcriptomic', 'radiomic', 'combined']):
    sv, X, tab, n, nb = out[key]; t = tab.head(15)[::-1]
    ax = fig.add_subplot(gs[row, 0])
    ax.barh(t.index, t['mean_abs_shap'], color=BLUE, height=0.62)
    for yi, (g, rr) in enumerate(t.iterrows()):
        ax.text(rr['mean_abs_shap'], yi, f"  {rr['pct_bootstrap_top10']:.0f}%", va='center', fontsize=7.5, color=INK2)
    ax.set_xlim(0, t['mean_abs_shap'].max() * 1.25); ax.set_xlabel('Mean |SHAP value|', fontsize=9)
    ax.tick_params(axis='y', length=0, labelsize=8); ax.spines['left'].set_visible(False)
    ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title(f'{titles[key]} (n = {n}, {nb} BCR)', fontsize=9, loc='left', color=INK)
    ax.text(-0.02, 1.05, 'ACE'[row], transform=ax.transAxes, fontsize=13, fontweight='bold', ha='right')
    ax2 = fig.add_subplot(gs[row, 1]); order = list(tab.index[:15]); col = [list(X.columns).index(g) for g in order]
    shap.plots.beeswarm(shap.Explanation(values=sv[:, col], data=X[order].values, feature_names=order),
                        max_display=15, order=list(range(15)), ax=ax2, show=False, plot_size=None)
    ax2.tick_params(axis='y', labelsize=8); ax2.tick_params(axis='x', labelsize=8)
    ax2.set_xlabel('SHAP value (impact on predicted BCR)', fontsize=9)
    ax2.text(-0.02, 1.05, 'BDF'[row], transform=ax2.transAxes, fontsize=13, fontweight='bold', ha='right')
fig.text(0.01, 0.005, 'Labels in A, C, E: % of 200 bootstrap refits in which the feature ranked among the 10 most important.',
         fontsize=8, color=INK2)
for ext in ('pdf', 'png'):
    fig.savefig(f'{FIG}/Figure4.{ext}', dpi=300, bbox_inches='tight')
print('Figure 4 written')
