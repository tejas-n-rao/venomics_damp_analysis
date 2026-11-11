"""Descriptive tolerances and measured within-species contrasts; no error attribution."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis_outputs'
COLS = ['apex_1_1_E._coli_ATCC_11775',
        'mbc_attention_E._coli_unspecified_strain',
        'ampredictor_E._coli_unspecified_strain']
raw = pd.read_csv(ROOT / 'df_damp_195625x57.tsv', sep='\t', usecols=COLS)
valid = np.isfinite(raw).all(axis=1) & raw.gt(0).all(axis=1)
z = np.log2(raw.loc[valid])
spread = z.max(axis=1) - z.min(axis=1)
# This cohort-specific transformation removes offsets, not experimental errors.
centered = z - z.median()
centered_spread = centered.max(axis=1) - centered.min(axis=1)
tolerances = pd.DataFrame([
    dict(maximum_fold_spread=t, raw_fraction=float(spread.le(np.log2(t)).mean()),
         centered_fraction=float(centered_spread.le(np.log2(t)).mean()), n=len(z))
    for t in [2, 4, 8, 16, 32, 64, 128, 256]
])
tolerances.to_csv(OUT / 'mic_variability_tolerance.csv', index=False)
summary = {'n': len(z), 'excluded': int((~valid).sum()),
           'median_raw_fold_spread': float(np.exp2(spread.median())),
           'median_centered_fold_spread': float(np.exp2(centered_spread.median())),
           'centering_medians_log2_uM': z.median().to_dict()}
(OUT / 'mic_variability_summary.json').write_text(json.dumps(summary, indent=2) + '\n')

# APEXDB has measured labels, unlike the prediction TSV. Keep only paired numeric
# positive labels; "inactive" has no verified threshold here and stays excluded.
db = pd.read_csv(ROOT / 'scratchpad/replay/APEXDB.csv')
assert not db.Sequence.duplicated().any()
pairs = [
    ('E. coli ATCC 11775', 'E. coli AIC221'),
    ('E. coli ATCC 11775', 'E. coli AIC222 - CRE'),
    ('E. coli AIC221', 'E. coli AIC222 - CRE'),
    ('P. aeruginosa PAO1', 'P. aeruginosa PA14'),
    ('S. aureus ATCC 12600', 'S. aureus ATCC BAA-1556 - MRSA'),
]
records = []
for a, b in pairs:
    x = db[[a, b]].apply(pd.to_numeric, errors='coerce')
    ok = np.isfinite(x).all(axis=1) & x.gt(0).all(axis=1)
    delta = np.log2(x.loc[ok, b]) - np.log2(x.loc[ok, a])
    absolute = delta.abs()
    records.append(dict(strain_a=a, strain_b=b, n=len(delta),
        excluded_from_1736=int((~ok).sum()),
        median_absolute_fold_ratio=float(np.exp2(absolute.median())),
        p90_absolute_fold_ratio=float(np.exp2(absolute.quantile(.9))),
        fraction_within_twofold=float(absolute.le(1).mean()),
        fraction_over_fourfold=float(absolute.gt(2).mean()),
        median_signed_log2_b_minus_a=float(delta.median())))
strain = pd.DataFrame(records)
strain.to_csv(OUT / 'mic_variability_strain_labels.csv', index=False)

fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
grid = np.linspace(0, 10, 401)
for values, label, color in [(spread, 'Exported predictions', '#b54a34'),
                            (centered_spread, 'After model-wise median centering', '#246b8e')]:
    sorted_values = np.sort(values)
    cdf = np.searchsorted(sorted_values, grid, side='right') / len(values)
    ax.plot(grid, 100 * cdf, label=label, color=color, lw=2)
ax.set(xticks=range(11), xticklabels=[str(2**i) for i in range(11)],
       xlabel='Allowed maximum/minimum prediction ratio (fold; log scale)',
       ylabel='Sequences within tolerance (%)', ylim=(0, 101), xlim=(0, 10),
       title=f'Three-family E. coli prediction spread (n = {len(z):,})')
ax.grid(alpha=.2)
ax.legend(loc='lower right')
fig.savefig(OUT / 'mic_variability_tolerance.png', dpi=180)
plt.close(fig)
print(json.dumps(summary, indent=2))
print(tolerances.to_string(index=False))
print(strain.to_string(index=False))
