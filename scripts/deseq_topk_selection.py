"""DESeq2-based gene pre-selection (top 10 / 20 / 30 by adjusted p-value) for BCR prediction.

Three settings, same classifiers and CV splits:
  inside   - DESeq2 re-run within each training fold; test patients never used for selection (reportable)
  outside  - DESeq2 on all 35 patients, then cross-validation (leakage; what "pre-selection" usually means)
  outside_random - as 'outside' but with permuted BCR labels (shows how much AUC leakage alone produces)
Classifiers: XGBoost with the settings used for the explanatory models, and ridge logistic regression (C = 0.1).
Usage: python deseq_topk_selection.py <out_dir>
"""
import sys, json, warnings
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')
from paths import DATA
OUT = sys.argv[1]; BASE = DATA
KS, REPS, N_PERM = (10, 20, 30), 5, 3

cnt = pd.read_csv(BASE + 'genes_with_ensembl_ids.csv', low_memory=False)
cnt = cnt.drop(columns=[c for c in cnt.columns if c in ('Ensembl ID', 'Unnamed: 36')]).drop_duplicates('Gene name').set_index('Gene name').apply(pd.to_numeric, errors='coerce').fillna(0)  # one copy per gene: BioMart merge duplicated rows
cnt.columns = [c.upper() for c in cnt.columns]
cnt = cnt[((cnt >= 10).mean(axis=1) >= 0.10)].astype(int)                 # label-free filter
lg = np.log(cnt.replace(0, np.nan)); ok = lg.notna().all(axis=1)
E = np.log2(cnt / np.exp(lg[ok].sub(lg[ok].mean(axis=1), axis=0).median()) + 1).T   # label-free normalisation
lab = pd.read_csv(BASE + 'radio_T2.csv', decimal=','); lab.index = lab['Samples'].str.replace('-', '').str.upper()
y_true = lab['Biochemical Recurrence'].reindex(E.index).astype(int).values
samples = E.index.values

def deseq_rank(idx, y):
    """Genes ranked by DESeq2 adjusted p-value (ties by |stat|) using only samples idx."""
    meta = pd.DataFrame({'BCR': y[idx].astype(str)}, index=samples[idx])
    dds = DeseqDataSet(counts=cnt[samples[idx]].T, metadata=meta, design='~BCR', quiet=True,
                       inference=DefaultInference(n_cpus=4))
    dds.deseq2()
    st = DeseqStats(dds, contrast=['BCR', '1', '0'], quiet=True, inference=DefaultInference(n_cpus=4)); st.summary()
    r = st.results_df.dropna(subset=['padj']).copy(); r['abs'] = r['stat'].abs()
    return r.sort_values(['padj', 'abs'], ascending=[True, False]).index.tolist()

CLFS = {'xgboost': lambda: XGBClassifier(n_estimators=200, max_depth=5, learning_rate=0.1, subsample=0.8,
                                          colsample_bytree=0.8, eval_metric='logloss', random_state=42, n_jobs=4),
        'ridge_logistic': lambda: LogisticRegression(C=0.1, max_iter=5000)}

def fold_predictions(y, tr, te, ranked):
    out = {}
    for k in KS:
        g = ranked[:k]; sc = StandardScaler().fit(E.iloc[tr][g])
        for name, mk in CLFS.items():
            m = mk().fit(sc.transform(E.iloc[tr][g]), y[tr])
            out[(k, name)] = m.predict_proba(sc.transform(E.iloc[te][g]))[:, 1]
    return out

def run_cv(y, rep, ranked_global=None):
    splits = list(StratifiedKFold(5, shuffle=True, random_state=rep).split(E, y))
    ranks = Parallel(n_jobs=5)(delayed(deseq_rank)(tr, y) for tr, _ in splits) if ranked_global is None else [ranked_global] * 5
    oof = {(k, n): np.zeros(len(y)) for k in KS for n in CLFS}
    for (tr, te), ranked in zip(splits, ranks):
        for key, p in fold_predictions(y, tr, te, ranked).items(): oof[key][te] = p
    return {key: roc_auc_score(y, p) for key, p in oof.items()}, [r[:30] for r in ranks]

rows, fold_genes = [], []
all_idx = np.arange(len(samples))
ranked_all = deseq_rank(all_idx, y_true)
for rep in range(REPS):
    auc_in, genes_in = run_cv(y_true, rep); fold_genes += genes_in
    auc_out, _ = run_cv(y_true, rep, ranked_all)
    for (k, n), a in auc_in.items(): rows.append(('inside', rep, k, n, a))
    for (k, n), a in auc_out.items(): rows.append(('outside', rep, k, n, a))
rng = np.random.default_rng(7)
for p in range(N_PERM):
    y_perm = rng.permutation(y_true)
    ranked_perm = deseq_rank(all_idx, y_perm)
    for rep in range(REPS):
        auc_p, _ = run_cv(y_perm, rep, ranked_perm)
        for (k, n), a in auc_p.items(): rows.append(('outside_random', p * REPS + rep, k, n, a))

df = pd.DataFrame(rows, columns=['setting', 'rep', 'top_k', 'classifier', 'auc'])
df.to_csv(f'{OUT}/deseq_topk_auc_all.csv', index=False)
summ = df.groupby(['setting', 'classifier', 'top_k'])['auc'].agg(
    mean='mean', p2_5=lambda s: np.percentile(s, 2.5), p97_5=lambda s: np.percentile(s, 97.5)).round(2)
summ.to_csv(f'{OUT}/deseq_topk_auc_summary.csv')
# stability: how often each gene of the full-cohort top 30 is re-selected inside the training folds
freq = pd.Series([g for gl in fold_genes for g in gl]).value_counts() / len(fold_genes)
stab = pd.DataFrame({'full_cohort_rank': range(1, 31)}, index=ranked_all[:30]); stab['selected_in_folds_top30'] = freq.reindex(stab.index).fillna(0).round(2)
stab.to_csv(f'{OUT}/deseq_topk_gene_stability.csv')
print(summ.to_string()); print(); print(stab.head(15).to_string())
