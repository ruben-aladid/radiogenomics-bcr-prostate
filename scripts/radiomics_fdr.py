# Table S2 (Mann-Whitney) and Table 2 (Fisher z) from radio_T2v1.csv. Usage: python radiomics_fdr.py <out_dir>
# and adds Benjamini-Hochberg FDR-adjusted p-values per family of tests.
import sys, pandas as pd, numpy as np, itertools
from scipy.stats import mannwhitneyu, spearmanr, norm
from statsmodels.stats.multitest import multipletests
from paths import DATA
base=DATA
feats=['SurfaceArea','HighGrayLevelRunEmphasis','Kurtosis','RobustMeanAbsoluteDeviation','LongRunHighGrayLevelEmphasis','ShortRunHighGrayLevelEmphasis','StandardDeviation','RunEntropy','Range','Maximum3DDiameter','MinorAxis','Maximum','MajorAxis','Autocorrelation','LongRunEmphasis','RunLengthNonUniformity','Maximum2DDiameterSlice','Energy','SumAverage','AverageIntensity','LeastAxis','GrayLevelNonUniformity','Minimum','Volume','RootMeanSquared']
d=pd.read_csv(base+'radio_T2v1.csv', decimal=',')
for c in d.columns[1:]: d[c]=pd.to_numeric(d[c].astype(str).str.replace(',','.'),errors='coerce')
g1=d[d['Biochemical Recurrence']==1]; g0=d[d['Biochemical Recurrence']==0]
mw=pd.DataFrame([(c,*mannwhitneyu(g1[c].dropna(),g0[c].dropna())) for c in feats],columns=['Variable','MannWhitneyU','p'])
mw['p_BH']=multipletests(mw.p,method='fdr_bh')[1]
mw.sort_values('p').to_csv(sys.argv[1]+'/radiomics_mannwhitney_fdr.csv',index=False)
rows=[]
for a,b in itertools.combinations(feats,2):
    x0=g0[[a,b]].dropna(); x1=g1[[a,b]].dropna()
    r0=spearmanr(x0[a],x0[b])[0]; r1=spearmanr(x1[a],x1[b])[0]
    z=(np.arctanh(r1)-np.arctanh(r0))/np.sqrt(1/(len(x1)-3)+1/(len(x0)-3))
    rows.append((a,b,r0,r1,z,2*(1-norm.cdf(abs(z)))))
F=pd.DataFrame(rows,columns=['Variable1','Variable2','rho_NED','rho_BCR','z','p']).dropna()
F['p_BH']=multipletests(F.p,method='fdr_bh')[1]
F.sort_values('p').to_csv(sys.argv[1]+'/radiomics_fisher_z_fdr.csv',index=False)
print(mw.sort_values('p').head(3).round(3).to_string()); print('Fisher z pairs',len(F),'nominal<0.05',(F.p<0.05).sum(),'BH<0.05',(F.p_BH<0.05).sum(),'min q',round(F.p_BH.min(),3))
