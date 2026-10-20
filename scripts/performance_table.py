"""Table S5 (supplementary): cross-validated discrimination of transcriptomic, radiomic and combined models.

Patients with both RNA-seq and T2 radiomics (n = 34; 11 BCR). One set of splits for every model:
5-fold stratified cross-validation repeated 10 times. Inside each training fold only:
  - RNA: DESeq2 (verified labels) on the training samples -> 20 genes with smallest adjusted p-value;
          features = log2 DESeq2-normalised counts, standardised on the training fold.
  - Radiomics: all 94 T2 features, mean imputation and standardisation fitted on the training fold.
  - Early fusion: the 20 RNA genes + 94 radiomic features in one model.
  - Late fusion: mean of the predicted probabilities of the RNA and radiomic models.
Classifiers: ridge logistic regression (C = 0.1) and XGBoost (settings of the original RNA model).
Reported: mean AUC and 2.5th-97.5th percentile range over the 10 repetitions (out-of-fold predictions).
Usage: python performance_table.py <results_dir>
"""
import sys, warnings
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

warnings.filterwarnings('ignore')
from paths import DATA
OUT = sys.argv[1]; BASE = DATA
REPS, K = 10, 20

cnt = pd.read_csv(BASE + 'genes_with_ensembl_ids.csv', low_memory=False)
cnt = cnt.drop(columns=[c for c in cnt.columns if c in ('Ensembl ID', 'Unnamed: 36')]).drop_duplicates('Gene name') \
         .set_index('Gene name').apply(pd.to_numeric, errors='coerce').fillna(0)
cnt.columns = [c.upper() for c in cnt.columns]
cnt = cnt[((cnt >= 10).mean(axis=1) >= 0.10)].astype(int)
r = pd.read_csv(BASE + 'radio_T2v1.csv', decimal=',')
r.index = r['Samples'].str.upper().replace({'P1755': 'P17551'})
R = r.drop(columns=['Samples', 'Biochemical Recurrence']).apply(
    lambda s: pd.to_numeric(s.astype(str).str.replace(',', '.'), errors='coerce')).dropna(axis=1, how='all')
ids = [i for i in R.index if i in cnt.columns]
cnt, R = cnt[ids], R.loc[ids]
lg = np.log(cnt.replace(0, np.nan)); ok = lg.notna().all(axis=1)
E = np.log2(cnt / np.exp(lg[ok].sub(lg[ok].mean(axis=1), axis=0).median()) + 1).T
lab = pd.read_csv(BASE + 'radio_T2.csv', decimal=','); lab.index = lab['Samples'].str.replace('-', '').str.upper()
y = lab['Biochemical Recurrence'].reindex(ids).astype(int).values
print(f'patients {len(ids)}, BCR {y.sum()}, genes {cnt.shape[0]}, radiomic features {R.shape[1]}')

def top_genes(tr):
    meta = pd.DataFrame({'BCR': y[tr].astype(str)}, index=np.array(ids)[tr])
    dds = DeseqDataSet(counts=cnt.iloc[:, tr].T, metadata=meta, design='~BCR', quiet=True, inference=DefaultInference(n_cpus=4))
    dds.deseq2()
    st = DeseqStats(dds, contrast=['BCR', '1', '0'], quiet=True, inference=DefaultInference(n_cpus=4)); st.summary()
    res = st.results_df.dropna(subset=['padj']).copy(); res['abs'] = res['stat'].abs()
    return res.sort_values(['padj', 'abs'], ascending=[True, False]).index[:K].tolist()

splits = [(rep, tr, te) for rep in range(REPS) for tr, te in StratifiedKFold(5, shuffle=True, random_state=rep).split(E, y)]
genes = Parallel(n_jobs=8)(delayed(top_genes)(tr) for _, tr, _ in splits)

CLF = {'Ridge logistic': lambda: LogisticRegression(C=0.1, max_iter=5000),
       'XGBoost': lambda: XGBClassifier(n_estimators=200, max_depth=5, learning_rate=0.1, subsample=0.8,
                                         colsample_bytree=0.8, eval_metric='logloss', random_state=42, n_jobs=8)}
oof = {(m, c): np.zeros((REPS, len(y))) for m in ['RNA (top-20 DESeq2 genes)', 'Radiomics (94 T2 features)',
                                                   'Early fusion (RNA + radiomics)', 'Late fusion (mean probability)'] for c in CLF}
for (rep, tr, te), g in zip(splits, genes):
    Xr_tr, Xr_te = E.iloc[tr][g], E.iloc[te][g]
    s = StandardScaler().fit(Xr_tr); Xr_tr, Xr_te = s.transform(Xr_tr), s.transform(Xr_te)
    imp = SimpleImputer(strategy='mean').fit(R.iloc[tr]); Xi_tr, Xi_te = imp.transform(R.iloc[tr]), imp.transform(R.iloc[te])
    s2 = StandardScaler().fit(Xi_tr); Xi_tr, Xi_te = s2.transform(Xi_tr), s2.transform(Xi_te)
    for c, mk in CLF.items():
        p_rna = mk().fit(Xr_tr, y[tr]).predict_proba(Xr_te)[:, 1]
        p_rad = mk().fit(Xi_tr, y[tr]).predict_proba(Xi_te)[:, 1]
        p_early = mk().fit(np.c_[Xr_tr, Xi_tr], y[tr]).predict_proba(np.c_[Xr_te, Xi_te])[:, 1]
        oof[('RNA (top-20 DESeq2 genes)', c)][rep, te] = p_rna
        oof[('Radiomics (94 T2 features)', c)][rep, te] = p_rad
        oof[('Early fusion (RNA + radiomics)', c)][rep, te] = p_early
        oof[('Late fusion (mean probability)', c)][rep, te] = (p_rna + p_rad) / 2

rows = []
for (m, c), P in oof.items():
    a = np.array([roc_auc_score(y, P[i]) for i in range(REPS)])
    rows.append({'Model': m, 'Classifier': c, 'AUC mean': round(a.mean(), 2),
                 'AUC 2.5%': round(np.percentile(a, 2.5), 2), 'AUC 97.5%': round(np.percentile(a, 97.5), 2)})
tab = pd.DataFrame(rows); tab.to_csv(f'{OUT}/table3_cv_performance.csv', index=False)
freq = pd.Series([x for g in genes for x in g]).value_counts() / len(genes)
freq.rename('pct_training_folds_top20').to_csv(f'{OUT}/table3_gene_selection_frequency.csv')
print(tab.to_string(index=False)); print(freq.head(10).round(2).to_string())
