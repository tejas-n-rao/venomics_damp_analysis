"""Exploratory near-neighbor audit; run from the project root."""
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
out=Path('analysis_outputs')
data=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
sources=pd.read_csv('scratchpad/source_sequences.tsv',sep='\t').set_index('sequence').source
cols=['apex_1_0_E._coli_ATCC11775','apex_1_1_E._coli_ATCC_11775','mbc_attention_E._coli_unspecified_strain','ampredictor_E._coli_unspecified_strain']
values=data[cols].dropna(); logs=np.log2(values)
spread=logs.max(axis=1)-logs.min(axis=1)
# Match each extreme sequence to a non-extreme sequence of the same source and length.
top=spread.nlargest(200).index.tolist(); top_set=set(top)
rng=np.random.default_rng(20261022)
pools=defaultdict(list)
for s in values.index:
    if s not in top_set: pools[(sources[s],len(s))].append(s)
controls=[]; match_rows=[]
for s in top:
    pool=pools[(sources[s],len(s))]
    if pool:
        control=pool.pop(int(rng.integers(len(pool))));controls.append(control)
        match_rows.append({'extreme':s,'control':control,'source':sources[s],'length':len(s)})
queries=top+controls
# Same length + same deletion position + same remaining sequence identifies
# one-substitution neighbors exactly. Store signatures only for the queries.
signatures=defaultdict(list)
for s in queries:
    for pos in range(len(s)):
        signatures[(len(s),pos,s[:pos]+s[pos+1:])].append(s)
neighbors={s:set() for s in queries}
for candidate in data.index:
    for pos in range(len(candidate)):
        for s in signatures.get((len(candidate),pos,candidate[:pos]+candidate[pos+1:]),()):
            if candidate!=s: neighbors[s].add(candidate)
rows=[];pairs=[]
for s in queries:
    rows.append({'sequence':s,'group':'top200' if s in top_set else 'matched_control',
        'source':sources[s],'length':len(s),'log2_disagreement':spread[s],
        'one_substitution_neighbors':len(neighbors[s]),
        'neighbors_in_top200':len(neighbors[s]&top_set)})
    if s in top_set:
        for t in sorted(neighbors[s]&top_set):
            if s<t: pairs.append({'sequence_a':s,'sequence_b':t,'length':len(s),
                'identity':(len(s)-1)/len(s),'range_a':spread[s],'range_b':spread[t]})
rows=pd.DataFrame(rows);rows.to_csv(out/'sequence_neighbor_audit.csv',index=False)
pd.DataFrame(match_rows).to_csv(out/'sequence_neighbor_controls.csv',index=False)
pd.DataFrame(pairs).to_csv(out/'extreme_one_substitution_pairs.csv',index=False)
summary=rows.groupby('group').agg(n=('sequence','size'),with_neighbor=('one_substitution_neighbors',lambda x:(x>0).sum()),median_neighbors=('one_substitution_neighbors','median'),with_extreme_neighbor=('neighbors_in_top200',lambda x:(x>0).sum()))
summary.to_csv(out/'sequence_neighbor_summary.csv');print(summary.to_string())
print('Pairs among the top 200:',len(pairs))
# Highlight examples whose relationships are directly inspectable; no motif annotation inferred.
examples=['VGIGTPISFYGGGAGHVPEYF','GGAGHVPEYFVGIGTPISFYG','GGAGHVPEYFVYIGTPISFYG','TTKNYGNGVCNSVNWCQCGNVWASCNLATGCAAWLCKLA','KYYGNGVCNSVNWCQCGNVWASCNLATGCAAWLCKLA']
e=values.loc[examples].copy(); e['source']=sources.loc[examples];e['log2_range']=spread.loc[examples]
e.to_csv(out/'sequence_similarity_examples.csv');print(e.to_string())
# Cheap global lower-bound audit: exact one-residue terminal extensions.
all_sequences=set(data.index)
terminal_pairs=[(s,s[1:]) for s in all_sequences if len(s)>=9 and s[1:] in all_sequences]
terminal_pairs += [(s,s[:-1]) for s in all_sequences if len(s)>=9 and s[:-1] in all_sequences]
terminal_pairs=sorted(set(terminal_pairs))
pd.DataFrame(terminal_pairs,columns=['longer','shorter']).to_csv(out/'terminal_extension_pairs.csv',index=False)
print('Unique one-residue terminal-extension pairs (shorter >=8 aa):',len(terminal_pairs))
print('Sequences participating:',len({s for pair in terminal_pairs for s in pair}))
