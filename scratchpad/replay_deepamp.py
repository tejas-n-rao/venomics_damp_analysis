from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,tensorflow as tf
O=Path('analysis_outputs');d=pd.read_csv('df_damp_195625x57.tsv',sep='\t').set_index('sequence')
# 32 seeded random eligible sequences plus 8 largest predictions, fixed before replay.
cols=[c for c in d if c.startswith('deep_amp')];ok=d[cols].dropna();s=list(dict.fromkeys(ok.sample(32,random_state=22).index.tolist()+ok.max(axis=1).nlargest(8).index.tolist()))
alphabet='BCDSQKIPTFNGHLRWAVEYM-';X=np.zeros((len(s),48,len(alphabet)),np.float32)
for i,seq in enumerate(s):
 for j,aa in enumerate(seq):X[i,j,alphabet.index(aa)]=1
T=X.copy()
for i,seq in enumerate(s):T[i,len(seq):,alphabet.index('-')]=1
rows=[]
for variant,prefix in [('CNN_gr_neg','deep_amp_cnn_neg'),('CNN_gr_pos','deep_amp_cnn_pos'),('LSTM_gr_neg','deep_amp_lstm_neg'),('LSTM_gr_pos','deep_amp_lstm_pos')]:
 m=tf.saved_model.load('scratchpad/replay/saved_models/paper/regressor/'+variant)
 fun=m.signatures['serving_default'];key=next(iter(fun.structured_input_signature[1]))
 col=next(c for c in cols if c.startswith(prefix))
 for mode,arr in [('zero_padding',X),('token_padding',T)]:
  pred=next(iter(fun(**{key:tf.constant(arr)}).values())).numpy().ravel()
  for seq,p in zip(s,pred): rows.append({'sequence':seq,'variant':variant,'padding':mode,'native_log10':float(p),'replayed_MIC':float(10.**p),'tsv_MIC':d.at[seq,col],'log2_ratio':float(p*np.log2(10)-np.log2(d.at[seq,col]))})
r=pd.DataFrame(rows);r.to_csv(O/'deepamp_replay.csv',index=False)
summary=r.groupby(['variant','padding']).log2_ratio.agg(n='size',median='median',max_abs=lambda x:x.abs().max(),median_abs=lambda x:x.abs().median());summary.to_csv(O/'deepamp_replay_summary.csv');print(summary.to_string())
files=[p for p in Path('scratchpad/replay/saved_models').rglob('*') if p.is_file()]
(O/'replay_manifest.json').write_text(json.dumps({'tensorflow':tf.__version__,'revision':'0a31ad796732b5e5874c3f92e771057e887f58a9','weights_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},'n_sequences':len(s)},indent=2))
