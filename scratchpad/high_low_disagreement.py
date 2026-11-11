"""Symmetric tails, matched contrasts, model direction and dependence sensitivities."""
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
O=Path('analysis_outputs');d=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
sources=pd.read_csv('scratchpad/source_sequences.tsv',sep='\t').set_index('sequence').source
cols={'APEX1.1':'apex_1_1_E._coli_ATCC_11775','MBC':'mbc_attention_E._coli_unspecified_strain','AMPredictor':'ampredictor_E._coli_unspecified_strain'}
v=d[list(cols.values())].rename(columns={b:a for a,b in cols.items()}).dropna();z=np.log2(v)
s=v.index.to_series();f=pd.DataFrame(index=s.index);f['length']=s.str.len()
f['charge_density']=(s.str.count('[KR]')-s.str.count('[DE]'))/f.length
f['hydrophobic_fraction']=s.str.count('[AILMFWVY]')/f.length;f['cys_fraction']=s.str.count('C')/f.length
f['entropy_bits']=[-sum((n/len(q))*np.log2(n/len(q)) for n in Counter(q).values()) for q in s]
f['source']=sources.reindex(f.index);f['log2_range']=z.max(axis=1)-z.min(axis=1);f['median_log2_MIC']=z.median(axis=1)
features=['length','charge_density','hydrophobic_fraction','cys_fraction','entropy_bits']
# Connected components of observed exact terminal extensions: a limited dependence sensitivity.
parent={q:q for q in d.index}
def root(q):
 while parent[q]!=q:parent[q]=parent[parent[q]];q=parent[q]
 return q
for q in d.index:
 if len(q)<9:continue
 for t in [q[1:],q[:-1]]:
  if t in parent:parent[root(q)]=root(t)
f['terminal_family']=[root(q) for q in f.index]
f['family_weight']=1/f.groupby('terminal_family').terminal_family.transform('size')
thresholds=[];summaries=[];matched=[];source_stats=[];directions=[];sensitivity=[];within=[]
labels={}
for pop,mask in [('all',pd.Series(True,index=f.index)),('restricted',f.length.between(8,39)&f.cys_fraction.eq(0))]:
 a=f.loc[mask].copy();lo,hi=a.log2_range.quantile([.05,.95]);a['group']=np.select([a.log2_range.le(lo),a.log2_range.ge(hi)],['low','high'],default='middle');labels[pop]=a.group
 thresholds.append({'population':pop,'n':len(a),'low_max_log2':lo,'high_min_log2':hi,'low_max_fold':2**lo,'high_min_fold':2**hi})
 for group,b in a.groupby('group'):
  rec={'population':pop,'group':group,'n':len(b),'median_log2_range':b.log2_range.median(),'median_MIC_uM':2**b.median_log2_MIC.median(),'cysteine_present_fraction':b.cys_fraction.gt(0).mean()}
  for feat in features:
   rec[feat+'_mean']=b[feat].mean();rec[feat+'_median']=b[feat].median();rec[feat+'_family_weighted_mean']=np.average(b[feat],weights=b.family_weight)
  summaries.append(rec)
  zz=z.loc[b.index];pair=zz.idxmax(axis=1)+' > '+zz.idxmin(axis=1)
  for name,count in pair.value_counts().items():directions.append({'population':pop,'group':group,'direction':name,'n':count,'fraction':count/len(b)})
  for source,g in b.groupby('source'):
   source_stats.append({'population':pop,'group':group,'source':source,'n':len(g),'fraction_of_group':len(g)/len(b),'fraction_of_source':len(g)/a.source.eq(source).sum()})
 # High versus remaining sequences, and direct high-versus-low controls.
 for contrast,control in [('high_vs_rest',a[a.group!='high']),('high_vs_low',a[a.group=='low'])]:
  pools=defaultdict(list);rng=np.random.default_rng(20261023)
  for q,r in control.iterrows():pools[(r.source,r.length,int(np.floor(r.median_log2_MIC)))].append(q)
  pairs=[]
  for q in rng.permutation(a.index[a.group.eq('high')]):
   r=a.loc[q];pool=pools[(r.source,r.length,int(np.floor(r.median_log2_MIC)))]
   if pool:pairs.append((q,pool.pop(int(rng.integers(len(pool))))))
  pd.DataFrame(pairs,columns=['high_sequence','control_sequence']).to_csv(O/f'{pop}_{contrast}_matches.csv',index=False)
  for feat in features:
   delta=np.array([a.at[q,feat]-a.at[t,feat] for q,t in pairs]);matched.append({'population':pop,'contrast':contrast,'feature':feat,'n_pairs':len(pairs),'high_n':a.group.eq('high').sum(),'mean_paired_difference':delta.mean() if len(delta) else np.nan,'median_paired_difference':np.median(delta) if len(delta) else np.nan})
 # Rank spread is invariant to monotonic rescaling; centering removes global offsets only.
 zz=z.loc[a.index];centered=zz-zz.median();cr=centered.max(axis=1)-centered.min(axis=1)
 ranks=zz.rank(pct=True);rr=ranks.max(axis=1)-ranks.min(axis=1)
 raw_high=set(a.index[a.group.eq('high')]);raw_low=set(a.index[a.group.eq('low')])
 for name,score in [('offset_centered',cr),('percentile_rank_range',rr)]:
  h=set(score.index[score.ge(score.quantile(.95))]);l=set(score.index[score.le(score.quantile(.05))])
  sensitivity.append({'population':pop,'definition':name,'high_retention':len(h&raw_high)/len(raw_high),'low_retention':len(l&raw_low)/len(raw_low),'spearman_with_raw':score.corr(a.log2_range,method='spearman')})
 # Equal tail fractions within each collection remove unequal source representation.
 for source,b in a.groupby('source'):
  l,h=b.log2_range.quantile([.05,.95]);low=b[b.log2_range.le(l)];high=b[b.log2_range.ge(h)]
  for feat in features:within.append({'population':pop,'source':source,'feature':feat,'n_low':len(low),'n_high':len(high),'mean_high_minus_low':high[feat].mean()-low[feat].mean()})
 f.loc[a.index].assign(group=a.group).join(v).to_csv(O/f'high_low_sequences_{pop}.csv')
for name,obj in [('high_low_thresholds',thresholds),('high_low_profiles',summaries),('high_low_matched',matched),('high_low_sources',source_stats),('high_low_model_direction',directions),('high_low_definition_sensitivity',sensitivity),('high_low_within_source',within)]:pd.DataFrame(obj).to_csv(O/(name+'.csv'),index=False)
# Classifier agreement is an independent output type, not experimental validation.
cc=['macrel_general_AMP','ai4amp_general_AMP','amp_clip_general_AMP','amplify_general_AMP']
c=d.loc[f.index,cc].dropna();votes=c.ge(.5).sum(axis=1);classification=[]
for group in ['low','middle','high']:
 idx=labels['all'].index[labels['all'].eq(group)].intersection(c.index);vv=votes.loc[idx]
 classification.append({'group':group,'n':len(idx),'all_positive':vv.eq(4).mean(),'all_negative':vv.eq(0).mean(),'mixed':vv.between(1,3).mean()})
pd.DataFrame(classification).to_csv(O/'high_low_classifier_agreement.csv',index=False)
# Compact publication-style panels; units and definitions appear on every axis.
S=pd.DataFrame(summaries);P=S[S.population=='all'].set_index('group').loc[['low','middle','high']]
fig,axs=plt.subplots(2,3,figsize=(12,7));colors=['#4477AA','#BBBBBB','#CC6677']
for ax,feat in zip(axs.flat,features+['median_log2_MIC']):
 if feat=='median_log2_MIC':ys=P.median_MIC_uM;label='Median predicted MIC (µM)'
 else:ys=P[feat+('_median' if feat=='length' else '_mean')];label=feat.replace('_',' ')+(' (median)' if feat=='length' else ' (mean)')
 ax.bar(P.index,ys,color=colors);ax.set_ylabel(label)
fig.suptitle('Three E. coli model families: lower 5%, middle 90%, upper 5%');fig.tight_layout();fig.savefig(O/'high_low_profiles.png',dpi=160);plt.close()
fig,axs=plt.subplots(1,2,figsize=(12,5));ss=pd.DataFrame(source_stats);g=ss[(ss.population=='all')&ss.group.isin(['low','high'])].pivot(index='source',columns='group',values='fraction_of_source');g.plot.barh(ax=axs[0],color={'low':'#4477AA','high':'#CC6677'});axs[0].set_xlabel('Fraction of each source in a global tail')
dd=pd.DataFrame(directions);g=dd[(dd.population=='all')&dd.group.isin(['low','high'])].pivot(index='direction',columns='group',values='fraction').fillna(0);g.plot.barh(ax=axs[1],color={'low':'#4477AA','high':'#CC6677'});axs[1].set_xlabel('Fraction of each tail');fig.tight_layout();fig.savefig(O/'high_low_sources_drivers.png',dpi=160);plt.close()
print(pd.DataFrame(thresholds).to_string(index=False));print(P.to_string());print(pd.DataFrame(matched).to_string(index=False));print(pd.DataFrame(sensitivity).to_string(index=False));print(pd.DataFrame(classification).to_string(index=False))
