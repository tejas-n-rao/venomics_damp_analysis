from google.colab import drive
import numpy as np
import os
import pandas as pd
# 1. Mount Google Drive to the default mount point
drive.mount("/content/drive")
project_dir="/content/drive/MyDrive/AMR/DAMP"
os.chdir(project_dir)

list_output=[f'./output/{fi}' for fi in os.listdir('./output/')]
list_output

df_apex10_raw=pd.read_csv('./output/df_apex_1_0.tsv',sep='\t')
df_apex11_raw=pd.read_csv('./output/df_apex_1_1.tsv',sep='\t')

df_apex10_raw.shape , df_apex11_raw.shape

list_verify=[f'./data/{fi}' for fi in os.listdir('./data_damp/') if 'verify' in fi]
list_verify

df_published_venomics=pd.read_csv('./data/df_verify_venomics.tsv',sep='\t')
df_published_venomics.shape

col_apex1_inf=np.concatenate([np.arange(8),np.array([26,27,28])])+1
mat_published_venomics=df_published_venomics.iloc[:,np.arange(1,12)].values
mat_inference_venomics=df_apex10_raw.loc[df_apex10_raw.sequence.isin(df_published_venomics.Sequence.values),
                                         df_apex10_raw.columns[col_apex1_inf]].values
max_deviation_venomics=np.max(np.abs(mat_inference_venomics-mat_published_venomics))
max_deviation_venomics

df_published_archaeomics=pd.read_csv('./data/df_verify_archaeomics.tsv',sep='\t')
df_published_archaeomics.shape

col_=['Sequence']
col_.extend( list(df_published_archaeomics.columns[1:].values))
df_published_archaeomics.columns=col_

mat_published_archaeomics=df_published_archaeomics.iloc[:,1:-1].values
mat_inferred_archaeomics=df_apex11_raw.loc[df_apex11_raw.sequence.isin(df_published_archaeomics.Sequence.values),df_apex11_raw.columns[1:]].values
mat_published_archaeomics, mat_inferred_archaeomics

error_archaemoics=np.mean(np.abs(mat_inferred_archaeomics-mat_published_archaeomics))
error_archaemoics

df_published_dad=pd.read_csv('./data_damp/df_verify_dad.tsv',sep='\t')
df_published_ampdiff=pd.read_csv('./data_damp/df_verify_AMPDiff.tsv',sep='\t')


df_inference_dad=df_apex11_raw.loc[df_apex11_raw.sequence.isin(df_published_dad.sequence.values),df_apex11_raw.columns[0:]]
df_inference_ampdiff=df_apex11_raw.loc[df_apex11_raw.sequence.isin(df_published_ampdiff.sequence.values),df_apex11_raw.columns[0:]]

df_published_dad.head(4)

df_inf_dad=df_inference_dad.set_index('sequence',inplace=False)
df_inf_dad.loc[df_published_dad.sequence.values]

mat_published_dad=df_published_dad.iloc[:,1:-8].values
mat_inf_dad=df_inf_dad.loc[df_published_dad.sequence.values].values
np.mean(np.abs(mat_published_dad-mat_inf_dad))

df_inf_ampdiff=df_inference_ampdiff.set_index('sequence',inplace=False)
mat_inference_ampdiff=df_inf_ampdiff.loc[df_published_ampdiff.sequence.values].values
mat_published_ampdiff=df_published_ampdiff.iloc[:,1:-8].values
np.mean(np.abs(mat_inference_ampdiff-mat_published_ampdiff))

df_damp=pd.read_csv('./data_damp/df_damp_190973x4.tsv',sep='\t')
df_apex10=pd.merge(df_damp.loc[:,[ 'source', 'model','sequence']],df_apex10_raw, on='sequence', how='inner')
APEX11_TO_APEX10 = [7, 0, 4, 5, 6, 1, 2, 3, 26, 27, 28]
col_34_strains=df_apex10.columns[3:].values
col_11_strains=df_apex10.columns[3:].values[APEX11_TO_APEX10]

(df_apex10.loc[df_apex10.source=='Extinctomics',col_34_strains].median(axis=1)<80).mean() , (df_apex10.loc[df_apex10.source=='Venomics',col_11_strains].median(axis=1)<32).mean()

df_apex11=pd.merge(df_damp.loc[:,[ 'source', 'model','sequence']],df_apex11_raw, on='sequence', how='inner')
df_apex11.head(4)

col_11_strains_apex11=df_apex11.columns[3:].values
print((df_apex11.loc[df_apex11.source=='Archaeomics',col_11_strains_apex11].mean(axis=1)<100).mean())

df_generated=df_apex11.loc[df_apex11.source=='AMPDIFF_synth50k',col_11_strains_apex11]
df_generated.reset_index(drop=True,inplace=True)
df_train=df_apex11.loc[df_apex11.source=='DRAMP_APD_DBAASP',col_11_strains_apex11]
df_train.reset_index(drop=True,inplace=True)

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



data = pd.concat(
    [
        pd.DataFrame({"Group": name, "MIC": df.mean(axis=1).to_numpy()})
        for name, df in {
            "Generated": df_generated,
            "Train": df_train,
        }.items()
    ],
    ignore_index=True,
)

sns.set_theme(style="ticks", context="talk")
fig, ax = plt.subplots(figsize=(5, 7))

sns.violinplot(
    data=data, x="Group", y="MIC", hue="Group",
    order=["Generated", "Train"],
    palette={"Generated": "#D87379", "Train": "#A5A5A5"},
    legend=False, saturation=1, linewidth=0,
    inner="box", cut=2, bw_method="scott",
    density_norm="width",
    inner_kws={"box_width": 8, "whis_width": 1.8, "marker": "_"},
    ax=ax,
)

ax.set(xlabel="", ylabel="Average MIC (µM)")
sns.despine(ax=ax)
fig.tight_layout()
plt.show()

# Exact statistics of the values plotted.
display(data.groupby("Group", sort=False)["MIC"].describe().round(2))