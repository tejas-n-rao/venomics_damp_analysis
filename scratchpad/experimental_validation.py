from pathlib import Path
from itertools import combinations
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
O=Path('analysis_outputs');d=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
rows=[]
for name,sc in [('extinctome_assays_2.xlsx','Seq'),('extinctome_assays_3.xlsx','Sequence')]:
 a=pd.read_excel('scratchpad/DAMP/data/data_damp/raw/'+name,header=1,keep_default_na=False)
 for _,r in a.iterrows():
  s=r[sc].strip();val=pd.to_numeric(r['E. coli ATCC11775'],errors='coerce')
  rows.append({'sequence':s,'table':name,'raw_MIC':r['E. coli ATCC11775'],'MIC_uM':val,'in_TSV':s in d.index,'strain':'E. coli ATCC11775','unit':'uM','medium':'LB','exposure_h':20,'inoculum_cells_ml':2e6,'modification_status':'not independently resolved','doi':'10.1038/s41551-024-01201-x'})
a=pd.DataFrame(rows);a.to_csv(O/'experimental_label_audit.csv',index=False)
a=a[a.in_TSV].drop_duplicates('sequence').set_index('sequence')
cols=['apex_1_0_E._coli_ATCC11775','apex_1_1_E._coli_ATCC_11775','mbc_attention_E._coli_unspecified_strain','ampredictor_E._coli_unspecified_strain'];names=['APEX1.0','APEX1.1','MBC','AMPredictor']
v=d.loc[a.index,cols];v.columns=names;a=a.join(v)
# Families: transitive graph using >=80% global identity for equal-length pairs,
# or exact containment of at least 8 residues with >=80% length coverage.
seqs=list(a.index);parent={s:s for s in seqs}
def root(s):
 while parent[s]!=s:s=parent[s]
 return s
for s,t in combinations(seqs,2):
 identity=sum(x==y for x,y in zip(s,t))/len(s) if len(s)==len(t) else 0
 containment=(s in t or t in s) and min(len(s),len(t))>=8 and min(len(s),len(t))/max(len(s),len(t))>=.8
 if identity>=.8 or containment:parent[root(s)]=root(t)
a['family']=[root(s) for s in a.index]
# Numeric-only primary analysis. N.A. is preserved as unresolved, never imputed.
e=a.dropna(subset=['MIC_uM']+names).copy();e['log2_observed']=np.log2(e.MIC_uM)
e['disagreement']=np.log2(e[['APEX1.1','MBC','AMPredictor']]).max(axis=1)-np.log2(e[['APEX1.1','MBC','AMPredictor']]).min(axis=1)
metrics=[]
for m in names:
 e[m+'_error']=np.log2(e[m])-e.log2_observed
 # Leave-one-family-out intercept calibration; fit only on other families.
 calibrated=[]
 for s,r in e.iterrows():
  train=e[e.family!=r.family];offset=(train.log2_observed-np.log2(train[m])).median()
  calibrated.append(np.log2(r[m])+offset)
 e[m+'_calibrated_log2']=calibrated
 err=e[m+'_error'];ce=e[m+'_calibrated_log2']-e.log2_observed
 metrics.append({'model':m,'n':len(e),'families':e.family.nunique(),'median_signed_log2_error':err.median(),'MAE_log2':err.abs().mean(),'within_twofold':err.abs().le(1).mean(),'spearman':e[m].corr(e.MIC_uM,method='spearman'),'disagreement_vs_absolute_error_rho':e.disagreement.corr(err.abs(),method='spearman'),'family_heldout_offset_MAE':ce.abs().mean()})
e.to_csv(O/'experimental_ecoli_matched.csv');pd.DataFrame(metrics).to_csv(O/'experimental_metrics.csv',index=False)
# Audit workbook: retain chemical modification text; do not silently strip it to claim matches.
w=pd.read_excel('wet_lab_validated_amp_sources.xlsx',sheet_name='Assay records',header=3);w.to_csv(O/'workbook_assay_inventory.csv',index=False)
print('Assay audit:',len(a),'matched sequences;',len(e),'numeric E. coli observations;',e.family.nunique(),'families')
print(pd.DataFrame(metrics).to_string(index=False));print(e[['MIC_uM','disagreement']+names].to_string())
fig,ax=plt.subplots(1,2,figsize=(11,4))
for m in names:ax[0].scatter(e.log2_observed,np.log2(e[m]),label=m,alpha=.7)
ax[0].plot([0,9],[0,9],'k--');ax[0].set(xlabel='Observed log2 MIC (µM)',ylabel='Predicted log2 MIC (µM)',title=f'Numeric observations only (n={len(e)})');ax[0].legend(fontsize=8)
for m in names:ax[1].scatter(e.disagreement,e[m+'_error'].abs(),label=m)
ax[1].set(xlabel='Three-family log2 disagreement',ylabel='Absolute log2 error',title='Selected cohort; exploratory association')
fig.tight_layout();fig.savefig(O/'experimental_validation.png',dpi=160);plt.close()
