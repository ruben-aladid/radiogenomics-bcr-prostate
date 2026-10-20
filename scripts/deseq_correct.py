"""Differential expression BCR vs NED (DESeq2, clinical labels). Usage: python deseq_correct.py <out_csv>"""
import pandas as pd, numpy as np, warnings, sys
warnings.filterwarnings('ignore')
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from paths import DATA
base=DATA
c=pd.read_csv(base+'genes_with_ensembl_ids.csv',low_memory=False).drop_duplicates('Gene name').set_index('Ensembl ID')  # one copy per gene: BioMart merge duplicated rows
names=c['Gene name']; c=c.drop(columns=[x for x in c.columns if x in ('Gene name','Unnamed: 36')]).apply(pd.to_numeric,errors='coerce').fillna(0).astype(int)
c=c[~c.index.duplicated()]
r=pd.read_csv(base+'radio_T2.csv',decimal=','); r['id']=r['Samples'].str.replace('-','').str.upper()
truth=dict(zip(r['id'],r['Biochemical Recurrence'].astype(int)))
meta=pd.DataFrame({'BCR':[str(truth[s]) for s in c.columns]},index=c.columns)
# filter as described in Methods: remove genes with <10 reads in >90% of samples
keep=((c>=10).mean(axis=1)>=0.10)
cf=c[keep]
print('samples',c.shape[1],'BCR',(meta.BCR=='1').sum(),'genes kept',keep.sum(),'of',len(c))
dds=DeseqDataSet(counts=cf.T,metadata=meta,design_factors='BCR',quiet=True)
dds.deseq2()
st=DeseqStats(dds,contrast=['BCR','1','0'],quiet=True); st.summary()
res=st.results_df.copy(); res['Gene name']=names.reindex(res.index).values
res.to_csv(sys.argv[1])
sig=res[(res.padj<0.05)&(res.log2FoldChange.abs()>1)]
print('DEGs padj<0.05 & |LFC|>1:',len(sig),' up:',(sig.log2FoldChange>0).sum(),' down:',(sig.log2FoldChange<0).sum())
print('padj<0.1:',(res.padj<0.1).sum(),' min padj:',res.padj.min())
print(sig.sort_values('padj').head(25)[['Gene name','baseMean','log2FoldChange','padj']].round(4).to_string())
