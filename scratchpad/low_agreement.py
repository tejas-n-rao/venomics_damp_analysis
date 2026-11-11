from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
O=Path('analysis_outputs');O.mkdir(exist_ok=True)
d=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
src=pd.read_csv('scratchpad/source_sequences.tsv',sep='\t').set_index('sequence').source
mapping={'APEX1.0':'apex_1_0_E._coli_ATCC11775','APEX1.1':'apex_1_1_E._coli_ATCC_11775','MBC':'mbc_attention_E._coli_unspecified_strain','AMPredictor':'ampredictor_E._coli_unspecified_strain'}
v=d[list(mapping.values())].rename(columns={b:a for a,b in mapping.items()}).dropna();z=np.log2(v)
s=v.index.to_series();f=pd.DataFrame(index=v.index)
f['length']=s.str.len();f['charge_density']=(s.str.count('[KR]')-s.str.count('[DE]'))/f.length
f['hydrophobic_fraction']=s.str.count('[AILMFWVY]')/f.length;f['cys_fraction']=s.str.count('C')/f.length
f['entropy_bits']=[-sum((n/len(seq))*np.log2(n/len(seq)) for n in Counter(seq).values()) for seq in s]
f['source']=src.reindex(f.index)
# Restriction is informed by the published AMPredictor training set, not a proven applicability domain.
domain=f.length.between(8,39)&f.cys_fraction.eq(0)
summary=[];feature_rows=[];match_rows=[];low_sets={};source_rows=[]
sets={'four_outputs':list(z),'three_families':['APEX1.1','MBC','AMPredictor'],'MBC_AMPredictor':['MBC','AMPredictor'],'APEX_same_strain':['APEX1.0','APEX1.1']}
features=['length','charge_density','hydrophobic_fraction','cys_fraction','entropy_bits']
for population,mask in [('all',pd.Series(True,index=f.index)),('restricted',domain)]:
 for comparison,cols in sets.items():
    a=z.loc[mask,cols];r=a.max(axis=1)-a.min(axis=1);cut=r.quantile(.05);low=r.index[r.le(cut)];other=r.index[r.gt(cut)]
    low_sets[(population,comparison)]=set(low)
    med=a.median(axis=1)
    summary.append({'population':population,'comparison':comparison,'n':len(a),'low_n':len(low),'low_5pct_max_log2_range':cut,'median_log2_range':r.median(),'within_twofold_fraction':r.le(1).mean(),'low_median_MIC_uM':2**med.loc[low].median(),'other_median_MIC_uM':2**med.loc[other].median()})
    for group,idx in [('low',low),('other',other)]:
      for feature in features:
       feature_rows.append({'population':population,'comparison':comparison,'group':group,'feature':feature,'median':f.loc[idx,feature].median(),'mean':f.loc[idx,feature].mean()})
    if comparison=='three_families':
      for source in f.source.unique():
       ix=f.index[mask&f.source.eq(source)];source_rows.append({'population':population,'source':source,'n':len(ix),'low_fraction':len(set(ix)&set(low))/len(ix) if len(ix) else np.nan})
      # Paired descriptive contrasts, matching source, exact length and 1-log2 median-MIC bins.
      pools=defaultdict(list);rng=np.random.default_rng(20261022)
      for seq in other: pools[(f.at[seq,'source'],f.at[seq,'length'],int(np.floor(med[seq])))].append(seq)
      matched=[]
      for seq in rng.permutation(low):
       pool=pools[(f.at[seq,'source'],f.at[seq,'length'],int(np.floor(med[seq])))]
       if pool: matched.append((seq,pool.pop(int(rng.integers(len(pool))))))
      for feat in features:
       delta=[f.at[x,feat]-f.at[y,feat] for x,y in matched]
       match_rows.append({'population':population,'feature':feat,'matched_pairs':len(matched),'low_total':len(low),'mean_paired_difference':np.mean(delta),'median_paired_difference':np.median(delta)})
      pd.DataFrame(matched,columns=['low_sequence','control_sequence']).to_csv(O/f'low_disagreement_matches_{population}.csv',index=False)
      f.loc[low].join(v.loc[low]).assign(log2_range=r.loc[low]).to_csv(O/f'low_disagreement_sequences_{population}.csv')
summary=pd.DataFrame(summary);summary.to_csv(O/'low_disagreement_summary.csv',index=False)
fr=pd.DataFrame(feature_rows);fr.to_csv(O/'low_disagreement_features.csv',index=False)
mr=pd.DataFrame(match_rows);mr.to_csv(O/'low_disagreement_matched.csv',index=False)
pd.DataFrame(source_rows).to_csv(O/'low_disagreement_by_source.csv',index=False)
print(summary.to_string(index=False));print(mr.to_string(index=False))
print(fr[(fr.population=='all')&(fr.comparison=='three_families')].to_string(index=False))
# Re-estimate strongest model disagreement claims on the restricted population.
from itertools import combinations
rows=[]
for pop,mask in [('all',pd.Series(True,index=f.index)),('restricted',domain)]:
 a=z.loc[mask].sort_index();k=int(np.ceil(len(a)*.01))
 for l,r in combinations(a,2):
  delta=a[r]-a[l]
  rows.append({'population':pop,'left':l,'right':r,'n':len(a),'spearman':a[l].corr(a[r],method='spearman'),'median_log2_ratio':delta.median(),'within_twofold':delta.abs().le(1).mean(),'top1pct_overlap':len(set(a[l].nsmallest(k).index)&set(a[r].nsmallest(k).index))/k})
pd.DataFrame(rows).to_csv(O/'domain_restricted_agreement.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
fig,ax=plt.subplots(1,2,figsize=(11,4))
for pop,mask in [('all',pd.Series(True,index=f.index)),('restricted',domain)]:
 a=z.loc[mask,['APEX1.1','MBC','AMPredictor']];r=a.max(axis=1)-a.min(axis=1)
 ax[0].hist(r,bins=np.linspace(0,16,65),density=True,histtype='step',label=f'{pop} (n={len(r):,})',linewidth=1.8)
ax[0].set(xlabel='Log2 range across three E. coli families',ylabel='Density');ax[0].legend()
g=fr[(fr.population=='all')&(fr.comparison=='three_families')&(fr.feature.isin(['charge_density','hydrophobic_fraction','cys_fraction']))]
g.pivot(index='feature',columns='group',values='mean').plot.bar(ax=ax[1]);ax[1].set(ylabel='Mean residue fraction / charge proxy',xlabel='');ax[1].tick_params(axis='x',rotation=20)
fig.tight_layout();fig.savefig(O/'low_disagreement_profiles.png',dpi=160);plt.close()
fig,axs=plt.subplots(1,2,figsize=(11,4))
ver=pd.read_csv(O/'apex_versions.csv');ver=ver[~ver.alias_unverified]
axs[0].barh(ver.target.str.replace('vancomycin-resistant ','',regex=False),ver.median_fold_ratio);axs[0].set_xlabel('Median MIC v1.1 / v1.0')
r=pd.DataFrame(rows);g=r[r.left.eq('APEX1.1')].pivot(index='right',columns='population',values='top1pct_overlap');g.plot.bar(ax=axs[1]);axs[1].set(ylabel='Fraction of top 1% shortlist shared',xlabel='Comparator to APEX1.1');axs[1].tick_params(axis='x',rotation=0)
fig.tight_layout();fig.savefig(O/'audit_disagreement_summary.png',dpi=160);plt.close()
