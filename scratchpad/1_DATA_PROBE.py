from google.colab import drive
import os
import pandas as pd
# 1. Mount Google Drive to the default mount point
drive.mount("/content/drive")
project_dir="/content/drive/MyDrive/AMR/DAMP"
os.chdir(project_dir)
!ls

# Make a data directory to store the assets

import os
from datetime import datetime

base_dir = "data"

if not os.path.exists(base_dir):
    os.makedirs(base_dir)
    target_dir = base_dir
else:
    date_str = datetime.now().strftime("%Y-%m-%d")
    target_dir = f"{base_dir}_{date_str}"
    os.makedirs(target_dir, exist_ok=True)

print(f"Target directory ready: {target_dir}")

# Loading .xlsx files directly into a pandas DataFrame via URL
url_ext1='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-024-01201-x/MediaObjects/41551_2024_1201_MOESM3_ESM.xlsx'

url_ext2='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-024-01201-x/MediaObjects/41551_2024_1201_MOESM4_ESM.xlsx'


url_ext3='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-024-01201-x/MediaObjects/41551_2024_1201_MOESM5_ESM.xlsx'

df_ext1 = pd.read_excel(url_ext1,header=1)
df_ext2 = pd.read_excel(url_ext2,header=1)
df_ext3 = pd.read_excel(url_ext3,header=1)
print(df_ext1.shape)
print(df_ext2.shape)
print(df_ext3.shape)

df_ext1.head(2)

df_ext2.head(2)

df_ext3.head(2)

df_ext1.to_csv('./data/df_extinctome.tsv',sep='\t',index=False)
df_ext2.to_csv('./data/df_extinctome_2.tsv',sep='\t',index=False)
df_ext3.to_csv('./data/df_extinctome_3.tsv',sep='\t',index=False)

df_damp=pd.DataFrame(columns=['sequence','source','model','paper'])
df_damp['sequence']=df_ext1.Seq.values
df_damp['source']='Extinctomics'
df_damp['model']='APEX_1.0'
df_damp['paper']='https://doi.org/10.1038/s41551-024-01201-x'
df_damp.shape

df_damp.to_csv('./data/df_damp_ext.tsv',sep='\t',index=False)

url_ven1='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-60051-6/MediaObjects/41467_2025_60051_MOESM3_ESM.csv'
url_ven2='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-60051-6/MediaObjects/41467_2025_60051_MOESM4_ESM.csv'
# url_ven3='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-60051-6/MediaObjects/41467_2025_60051_MOESM7_ESM.xlsx' # 403 errors

url_ven3='https://github.com/vinayprabhu/DAMP/raw/refs/heads/main/data/41467_2025_60051_MOESM7_ESM.xlsx'

# Fix url_ven2 reading error by skipping malformed/unbalanced lines
df_ven1 = pd.read_csv(url_ven1)
# 1. Download raw CSV locally
!wget -q -O raw_ven2.csv "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-60051-6/MediaObjects/41467_2025_60051_MOESM4_ESM.csv"

# 2. Install csvkit for file cleanup
!pip install -q csvkit

# 3. Clean unbalanced quotes and delimiters without dropping any lines
!csvclean raw_ven2.csv

# 4. Read the cleaned output file
df_ven2 = pd.read_csv("raw_ven2.csv")
df_ven3 = pd.read_excel(url_ven3, engine="openpyxl")



print(df_ven1.shape)
print(df_ven2.shape)
print(df_ven3.shape)

df_ven2.head(2)

df_ven3.head(4)

df_ven1.to_csv('./data/df_verify_venomics.tsv',sep='\t',index=False)

df_damp_v=pd.DataFrame(columns=['sequence','source','model','paper'])
df_damp_v['sequence']=df_ven1.Sequence.values
df_damp_v['source']='Venomics'
df_damp_v['model']='APEX_1.0'
df_damp_v['paper']='https://doi.org/10.1038/s41467-025-60051-6'

df_damp_v.to_csv('./data/df_damp_v.tsv',sep='\t',index=False)

print(f'df_damp.shape before venomics merger: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_v])
print(f'df_damp.shape after venomics merger: {df_damp.shape}')
df_damp.shape

df_damp.paper.value_counts()

df_damp.source.value_counts()

url_arch1='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM3_ESM.csv'
url_arch2='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM4_ESM.csv'
url_arch3='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM5_ESM.csv'
url_arch4='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM6_ESM.csv'

import csv
# Specialized CSV parsing function after encountering a few errors
def read_csv(url):
  # Download the raw file
  !wget -q -O temp_csv.csv $url

  # Robustly parse CSV while retaining multi-line quoted fields
  rows = []
  with open("temp_csv.csv", mode="r", encoding="utf-8-sig", errors="replace") as f:
      reader = csv.reader(f)
      header = next(reader)
      for r in reader:
          rows.append(r)

  # Convert directly into pandas DataFrame
  df_out = pd.DataFrame(rows, columns=header)
  return df_out

df_arch1=read_csv(url_arch1)
df_arch2=read_csv(url_arch2)
df_arch3=read_csv(url_arch3)
df_arch4=read_csv(url_arch4)

df_arch1.shape,df_arch2.shape,df_arch3.shape,df_arch4.shape

df_arch1.head()

df_arch1.to_csv('./data/df_verify_archaeomics.tsv',sep='\t',index=False)

df_damp_arch=pd.DataFrame(columns=df_damp.columns)
df_damp_arch['sequence']=df_arch1.iloc[:,0].values
df_damp_arch['model']='APEX_1.1'
df_damp_arch['source']='Archaeomics'
df_damp_arch['paper']='https://doi.org/10.1038/s41564-025-02061-0'

df_damp_arch.head(4)

df_damp_arch.to_csv('./data/df_damp_arch.tsv',sep='\t',index=False)

print(f'Pre-Archaeomics merger shape: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_arch])
print(f'Post-Archaeomics merger shape: {df_damp.shape}')

df_damp.paper.value_counts()

df_damp.source.value_counts()

url_train_ampdiff_dad='https://raw.githubusercontent.com/vinayprabhu/DAMP/refs/heads/main/data/ampdiff_mmc2.csv'
url_ad='https://raw.githubusercontent.com/vinayprabhu/DAMP/refs/heads/main/data/ampdiff_mmc3.csv'

df_dad= pd.read_csv(url_train_ampdiff_dad)
df_dad.shape

df_dad.head(4)

df_damp_dad=pd.DataFrame(columns=df_damp.columns)
df_damp_dad['sequence']=df_dad.sequence.values
df_damp_dad['model']='CANONICAL_DATASETS'
df_damp_dad['source']='DRAMP_APD_DBAASP'
df_damp_dad['paper']='https://doi.org/10.1016/j.celbio.2025.100183'
print(f'Pre-merger shape: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_dad])
print(f'Post-merger shape: {df_damp.shape}')

df_damp_dad.to_csv('./data/df_damp_dad.tsv',sep='\t',index=False)

df_ampdiff=pd.read_csv(url_ad)
df_ampdiff.head(4)

df_ampdiff.to_csv('./data/df_verify_AMPDiff.tsv',sep='\t',index=False)

df_damp_ampdiff=pd.DataFrame(columns=df_damp.columns)
df_damp_ampdiff['sequence']=df_ampdiff.sequence.values
df_damp_ampdiff['model']='APEX_1.1'
df_damp_ampdiff['source']='AMPDIFF_synth50k'
df_damp_ampdiff['paper']='https://doi.org/10.1016/j.celbio.2025.100183'
df_damp_ampdiff.shape

df_damp_ampdiff.to_csv('./data/df_damp_ampdiff_50k.tsv',sep='\t',index=False)

print(f'Pre-merger shape: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_ampdiff])
print(f'Post-merger shape: {df_damp.shape}')

df_damp.source.value_counts()

!pip install -q biopython
import pandas as pd
from Bio import SeqIO

def fasta_to_df(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append({
            "ID": record.id,
            "Description": record.description,
            "Sequence": str(record.seq)
        })
    return pd.DataFrame(records)

# Use the FASTA filename specified in your active cell
file_human_proteome = 'uniprot_sprot_HomoSapiens_WithIsoF_2020_07_17_42361.fasta'
df_fasta = fasta_to_df(file_human_proteome)

print(df_fasta.shape)
df_fasta.head(3)

url_hp_1='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41551-021-00801-1/MediaObjects/41551_2021_801_MOESM3_ESM.csv'
df_hp=pd.read_csv(url_hp_1)
print(df_hp.shape)

df_hp.head()

import random
for trials in range(4):
  IND_TEST = random.randint(0, df_hp.shape[0])
  print(f'IND_TEST = {IND_TEST}')
  raw_seq=df_fasta.loc[df_fasta.ID == df_hp.loc[IND_TEST,'RecID'],'Sequence'].values[0]
  ind_beg=df_hp.loc[IND_TEST,'start']
  ind_end=ind_beg+df_hp.loc[IND_TEST,'len']
  if(raw_seq[ind_beg:ind_end]==df_hp.loc[IND_TEST,'sequence']):
    print(f'They match! {raw_seq[ind_beg:ind_end]}')
  else:
    print(f'{raw_seq[ind_beg:ind_end],df_hp.loc[IND_TEST,'sequence']}')


df_damp_hp=pd.DataFrame(columns=df_damp.columns)
df_damp_hp['sequence']=df_hp.sequence.values
df_damp_hp['model']='AMP_potency_score_pane_et_al'
df_damp_hp['source']='Human_Proteomics'
df_damp_hp['paper']='https://doi.org/10.1038/s41551-021-00801-1'
df_damp_hp.shape

df_damp_hp.to_csv('./data/df_damp_hp.tsv',sep='\t',index=False)
print(f'Pre-merger shape: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_hp])
print(f'Post-merger shape: {df_damp.shape}')

url_fasta= 'https://raw.githubusercontent.com/szczurek-lab/amp-challenge-2027/refs/heads/main/data/antibacterial.fasta'

import requests
import io
import pandas as pd
from Bio import SeqIO

# Fetch the fasta file from the URL
response = requests.get(url_fasta)
response.raise_for_status()

# Parse the text stream and extract records
records = []
for record in SeqIO.parse(io.StringIO(response.text), "fasta"):
    records.append({
        "ID": record.id,
        "Description": record.description,
        "Sequence": str(record.seq)
    })

# Convert directly into a pandas DataFrame
df_neurips = pd.DataFrame(records)
print(df_neurips.shape)
df_neurips.head(2)

df_damp_neurips=pd.DataFrame(columns=df_damp.columns)
df_damp_neurips['sequence']=df_neurips.Sequence.values
df_damp_neurips['model']='UNCLEAR'
df_damp_neurips['source']='NEURIPS_AMP_CHALLENGE'
df_damp_neurips['paper']='https://github.com/szczurek-lab/amp-challenge-2027/blob/main/AMPChallenge_NeurIPSCompetition.pdf'
df_damp_neurips.shape

df_damp_neurips.to_csv('./data/df_damp_neurips.tsv',sep='\t',index=False)
print(f'Pre-merger shape: {df_damp.shape}')
df_damp=pd.concat([df_damp,df_damp_neurips])
print(f'Post-merger shape: {df_damp.shape}')

print(df_damp.shape)
df_damp.to_csv('./data/df_damp_09_16_26.tsv',sep='\t',index=False)
df_damp.source.value_counts()

df_damp.model.value_counts()

(df_damp.sequence.value_counts()>1).sum(),(df_damp.sequence.value_counts()>2).sum(),df_damp.shape[0]-df_damp.sequence.nunique()

df_damp = df_damp.drop_duplicates(subset=['sequence'], keep='first')
# Reset the index and drop the old index column
df_damp = df_damp.reset_index(drop=True)
df_damp.source.value_counts()

df_damp.to_csv('./data/df_damp_190973x4.tsv',sep='\t', index=False)

df_d2=pd.read_csv('./data_damp/df_damp_190973x4.tsv',sep='\t')
df_d2.shape

df_d2.equals(df_damp)

df_d1=pd.read_csv('./data/df_damp_09_16_26.tsv',sep='\t')
df_d2=pd.read_csv('./data/df_damp_190973x4.tsv',sep='\t')
df_d1.shape,df_d2.shape

df_d1

# 1. Identify sequences appearing more than once
counts = df_d1['sequence'].value_counts()
dupes = counts[counts > 1].index

# 2. Filter the original dataframe (using df_d1 based on your active cell)
df_filtered = df_d1[df_d1['sequence'].isin(dupes)]

# 3. Group by sequence and aggregate the source column with unique values
df_repeating = df_filtered.groupby('sequence', as_index=False).agg(
    source=('source', lambda x: ' | '.join(x.dropna().unique().astype(str)))
)

df_repeating.head()

df_repeating.source.value_counts()

df_damp.source.unique()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

# Map each source to its set of unique sequences using your base dataframe
source_sets = df_d1.groupby('source')['sequence'].apply(set).to_dict()
sources = sorted(source_sets.keys())

# Initialize the square overlap matrix
df_heatmap = pd.DataFrame(index=sources, columns=sources, dtype=float)

# Populate the matrix with set intersections and total lengths
for i, src1 in enumerate(sources):
    for j, src2 in enumerate(sources):
        if i == j:
            df_heatmap.loc[src1, src2] = len(source_sets[src1])
        else:
            overlap = len(source_sets[src1].intersection(source_sets[src2]))
            df_heatmap.loc[src1, src2] = overlap

# Generate a mask for the lower triangle
mask = np.tril(np.ones_like(df_heatmap, dtype=bool), k=-1)

# Create a plain white colormap
white_cmap = ListedColormap(['white'])

# Render the heatmap with a white background, black text, and cell borders
plt.figure(figsize=(12, 10))
sns.heatmap(
    df_heatmap,
    mask=mask,
    annot=True,
    fmt="g",
    cmap=white_cmap,
    cbar=False,
    linewidths=1,
    linecolor='black',
    annot_kws={"size": 10, "weight": "bold", "color": "black"}
)
plt.title("Total Sequences (Diagonal) and Overlaps (Off-Diagonal)")
plt.tight_layout()
plt.show()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import warnings
from matplotlib.colors import ListedColormap

# Suppress missing font warnings from matplotlib
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib.font_manager")

# Map each source to its set of unique sequences
source_sets = df_d1.groupby('source')['sequence'].apply(set).to_dict()
sources = sorted(source_sets.keys())

# Initialize the square overlap matrix
df_heatmap = pd.DataFrame(index=sources, columns=sources, dtype=float)

# Populate the matrix with set intersections and total lengths
for i, src1 in enumerate(sources):
    for j, src2 in enumerate(sources):
        if i == j:
            df_heatmap.loc[src1, src2] = len(source_sets[src1])
        else:
            overlap = len(source_sets[src1].intersection(source_sets[src2]))
            df_heatmap.loc[src1, src2] = overlap

# Generate a mask for the lower triangle
mask = np.tril(np.ones_like(df_heatmap, dtype=bool), k=-1)
white_cmap = ListedColormap(['white'])

# Apply XKCD styling context manager
with plt.xkcd():
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        df_heatmap,
        mask=mask,
        annot=True,
        fmt="g",
        cmap=white_cmap,
        cbar=False,
        linewidths=0.5,
        linecolor='black',
        annot_kws={"size": 11, "color": "black"}
    )
    plt.title("Total Sequences (Diagonal) and Overlaps (Off-Diagonal)")
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

