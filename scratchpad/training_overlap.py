from pathlib import Path
import numpy as np,pandas as pd
O=Path('analysis_outputs');d=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
t=pd.read_csv('scratchpad/replay/APEXDB.csv');t.Sequence=t.Sequence.str.strip()
e=pd.read_csv(O/'experimental_ecoli_matched.csv');r=e[['sequence','MIC_uM']].merge(t,left_on='sequence',right_on='Sequence');r[['sequence','MIC_uM','E. coli ATCC 11775','Group (APEX training)']].to_csv(O/'experimental_training_overlap.csv',index=False)
# Literal backslash-t separators in the pinned AMPredictor file.
lines=Path('scratchpad/replay/ampredictor_train.txt').read_text().splitlines();amp=set()
for line in lines:
 parts=line.replace('\\t','\t').split('\t')
 if len(parts)>=2:amp.add(parts[1].strip())
rows=[]
for name,seqs in [('APEX_CV',set(t.loc[t['Group (APEX training)'].eq('CV'),'Sequence'])),('APEX_Ind',set(t.loc[t['Group (APEX training)'].eq('Ind'),'Sequence'])),('AMPredictor_train',amp)]:
 rows.append({'corpus':name,'unique_corpus_sequences':len(seqs),'exact_TSV_overlap':len(set(d.index)&seqs),'numeric_assay_overlap':len(set(e.sequence)&seqs)})
pd.DataFrame(rows).to_csv(O/'training_overlap_counts.csv',index=False);print(pd.DataFrame(rows).to_string(index=False))
c=pd.read_csv(O/'extreme_sequence_components.csv');target=c.loc[c.component.eq(1),'sequence'];nearest=[]
for s in target:
 for name,seqs in [('APEX_corpus',set(t.Sequence)),('AMPredictor_train',amp)]:
  candidates=[q for q in seqs if len(q)==len(s)]
  if candidates:
   q=max(sorted(candidates),key=lambda q:sum(a==b for a,b in zip(s,q)))
   nearest.append({'sequence':s,'corpus':name,'nearest_equal_length_sequence':q,'ungapped_identity':sum(a==b for a,b in zip(s,q))/len(s),'has_cysteine':'C' in s})
pd.DataFrame(nearest).to_csv(O/'discordant_family_training_neighbors.csv',index=False)
print(pd.DataFrame(nearest).groupby('corpus').ungapped_identity.agg(['min','median','max']))
