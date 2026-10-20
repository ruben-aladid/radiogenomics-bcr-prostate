"""Descriptive (non-predictive) top-20 DESeq2 gene model with SHAP, verified BCR labels.

- Genes: 20 smallest DESeq2 adjusted p-values in the full cohort (ties by |stat|), one copy per gene.
- Model: XGBoost with the settings used for the explanatory models, fitted on all 35 patients; SHAP on all 35 patients.
  Because selection and fitting use every patient, no performance metric is derived from this model.
- Stability: DESeq2 re-run within the training folds of 5x repeated stratified 5-fold CV (25 runs);
  for each gene, the % of training folds in which it ranks in the top 20.
Usage: python top20_explanatory.py <results_dir>
"""
import sys, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shap
from joblib import Parallel, delayed
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from xgboost import XGBClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings('ignore')
from paths import DATA
RC = sys.argv[1]; EX = f'{RC}/explanatory'; BASE = DATA
cnt = pd.read_csv(BASE + 'genes_with_ensembl_ids.csv', low_memory=False)
cnt = cnt.drop(columns=[c for c in cnt.columns if c in ('Ensembl ID', 'Unnamed: 36')]).drop_duplicates('Gene name') \
         .set_index('Gene name').apply(pd.to_numeric, errors='coerce').fillna(0)
cnt.columns = [c.upper() for c in cnt.columns]
cnt = cnt[((cnt >= 10).mean(axis=1) >= 0.10)].astype(int)
lg = np.log(cnt.replace(0, np.nan)); ok = lg.notna().all(axis=1)
E = np.log2(cnt / np.exp(lg[ok].sub(lg[ok].mean(axis=1), axis=0).median()) + 1).T
lab = pd.read_csv(BASE + 'radio_T2.csv', decimal=','); lab.index = lab['Samples'].str.replace('-', '').str.upper()
y = lab['Biochemical Recurrence'].reindex(E.index).astype(int).values
samples = E.index.values

def deseq(idx):
    meta = pd.DataFrame({'BCR': y[idx].astype(str)}, index=samples[idx])
    dds = DeseqDataSet(counts=cnt[samples[idx]].T, metadata=meta, design='~BCR', quiet=True, inference=DefaultInference(n_cpus=4))
    dds.deseq2()
    st = DeseqStats(dds, contrast=['BCR', '1', '0'], quiet=True, inference=DefaultInference(n_cpus=4)); st.summary()
    r = st.results_df.dropna(subset=['padj']).copy(); r['abs'] = r['stat'].abs()
    return r.sort_values(['padj', 'abs'], ascending=[True, False])

full = deseq(np.arange(len(y)))
top = full.head(20)
splits = [tr for rep in range(5) for tr, _ in StratifiedKFold(5, shuffle=True, random_state=rep).split(E, y)]
ranks = Parallel(n_jobs=6)(delayed(lambda tr: deseq(tr).index[:20].tolist())(tr) for tr in splits)
freq = pd.Series([g for r in ranks for g in r]).value_counts() / len(ranks)

X = pd.DataFrame(StandardScaler().fit_transform(E[top.index]), index=E.index, columns=top.index)
m = XGBClassifier(n_estimators=200, max_depth=5, learning_rate=0.1, subsample=0.8, colsample_bytree=0.8,
                  eval_metric='logloss', random_state=42).fit(X, y)
sv = shap.TreeExplainer(m).shap_values(X)
tab = pd.DataFrame({'log2FC_BCR_vs_NED': top['log2FoldChange'].round(2), 'padj': top['padj'].round(4),
                    'mean_abs_shap': np.abs(sv).mean(axis=0).round(3),
                    'pct_training_folds_top20': (freq.reindex(top.index).fillna(0) * 100).round(0)})
tab = tab.sort_values('mean_abs_shap', ascending=False)
tab.to_csv(f'{EX}/top20_explanatory_table.csv')

INK, INK2, GRID, BLUE = '#0b0b0b', '#52514e', '#e6e5e0', '#2a78d6'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'pdf.fonttype': 42})
fig = plt.figure(figsize=(12, 6.2)); gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.25], wspace=0.55)
ax = fig.add_subplot(gs[0, 0]); t = tab[::-1]
ax.barh(t.index, t['mean_abs_shap'], color=BLUE, height=0.62)
for yi, (g, r) in enumerate(t.iterrows()):
    ax.text(r['mean_abs_shap'], yi, f"  {r['pct_training_folds_top20']:.0f}%", va='center', fontsize=7, color=INK2)
ax.set_xlim(0, t['mean_abs_shap'].max() * 1.25); ax.set_xlabel('Mean |SHAP value|')
ax.tick_params(axis='y', length=0, labelsize=8); ax.spines['left'].set_visible(False)
ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title('Top-20 DESeq2 genes (label: % of training folds in top 20)', fontsize=8.5, loc='left', color=INK)
ax.text(-0.02, 1.04, 'A', transform=ax.transAxes, fontsize=13, fontweight='bold', ha='right')
ax2 = fig.add_subplot(gs[0, 1])
order = list(tab.index)
shap.plots.beeswarm(shap.Explanation(values=sv[:, [list(X.columns).index(g) for g in order]], data=X[order].values,
                                     feature_names=order), max_display=20, order=list(range(20)), ax=ax2, show=False, plot_size=None)
ax2.tick_params(axis='y', labelsize=8); ax2.set_xlabel('SHAP value (impact on predicted BCR)')
ax2.text(-0.02, 1.04, 'B', transform=ax2.transAxes, fontsize=13, fontweight='bold', ha='right')
for ext in ('pdf', 'png'):
    fig.savefig(f'{EX}/Figure_top20_explanatory.{ext}', dpi=300, bbox_inches='tight')
print(tab.to_string())
