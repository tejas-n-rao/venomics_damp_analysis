"""Independent checks and focused extensions of DAMP's descriptive analysis.

Run from any directory with the existing analysis environment:
    MPLCONFIGDIR=/private/tmp/venomics-mpl python analysis_review/review_analysis.py

No network or model inference. Original inputs and outputs are read-only.
All new estimates are exploratory summaries of this finite, selected collection.
"""
from pathlib import Path
from collections import Counter
from itertools import combinations
import hashlib
import json
import os
import platform
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import scipy
from scipy.stats import kendalltau
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
COLS = {
    "APEX 1.0": "apex_1_0_E._coli_ATCC11775",
    "APEX 1.1": "apex_1_1_E._coli_ATCC_11775",
    "MBC-Attention": "mbc_attention_E._coli_unspecified_strain",
    "AMPredictor": "ampredictor_E._coli_unspecified_strain",
}
PRIMARY = ["APEX 1.1", "MBC-Attention", "AMPredictor"]
OLD_NAMES = {"APEX1.0": "APEX 1.0", "APEX1.1": "APEX 1.1",
             "MBC": "MBC-Attention", "AMPredictor": "AMPredictor"}
LONG_OLD = dict(zip(["APEX_1.0_ATCC11775", "APEX_1.1_ATCC11775", "MBC_Attention", "AMPredictor"], COLS))
checks = []


def check(name, observed, expected, rtol=1e-8, atol=1e-9):
    if isinstance(observed, (str, bool, np.bool_)):
        passed = bool(observed == expected)
    else:
        passed = bool(np.allclose(observed, expected, rtol=rtol, atol=atol, equal_nan=True))
    checks.append({"check": name, "passed": passed, "observed": str(observed), "expected": str(expected)})


def save(frame, name):
    frame.to_csv(OUT / f"{name}.csv", index=False)


def figsave(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=190, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def pair_metrics(z, population):
    rows = []
    for left, right in combinations(z, 2):
        delta = z[left] - z[right]
        residual = delta - delta.median()
        # Kendall tau-a: ties contribute zero; denominator is all possible pairs.
        # tau-b alone cannot be read directly as a pair-reversal proportion.
        tau_b = float(kendalltau(z[left], z[right]).statistic)
        n = len(z)
        pairs = n * (n - 1) / 2
        def tied_pairs(x):
            counts = x.value_counts().to_numpy(dtype=np.float64)
            return float(np.sum(counts * (counts - 1) / 2))
        ties_left, ties_right = tied_pairs(z[left]), tied_pairs(z[right])
        joint_counts = z.groupby([left, right]).size().to_numpy(dtype=np.float64)
        ties_joint = float(np.sum(joint_counts * (joint_counts - 1) / 2))
        comparable = pairs - ties_left - ties_right + ties_joint
        concordant_minus_discordant = tau_b * np.sqrt((pairs - ties_left) * (pairs - ties_right))
        discordant = (comparable - concordant_minus_discordant) / 2
        rows.append({"population": population, "left": left, "right": right, "n": n,
                     "spearman": z[left].corr(z[right], method="spearman"),
                     "median_log2_left_minus_right": delta.median(),
                     "median_left_right_fold_ratio": np.exp2(delta.median()),
                     "raw_median_absolute_log2_difference": delta.abs().median(),
                     "offset_removed_median_absolute_log2_difference": residual.abs().median(),
                     "raw_within_twofold": delta.abs().le(1).mean(),
                     "offset_removed_within_twofold": residual.abs().le(1).mean(),
                     "kendall_tau_b": tau_b, "strictly_comparable_pair_fraction": comparable / pairs,
                     "opposite_order_fraction_of_comparable_pairs": discordant / comparable})
    return rows


def main():
    path = ROOT / "df_damp_195625x57.tsv"
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    original_manifest = json.loads((ROOT / "analysis_outputs/manifest.json").read_text())
    check("input SHA256", checksum, original_manifest["input_sha256"])
    raw = pd.read_csv(path, sep="\t").set_index("sequence", verify_integrity=True)
    check("rows", len(raw), 195625)
    check("physical columns", len(raw.columns) + 1, 58)
    metadata = pd.DataFrame(json.loads((ROOT / "scratchpad/damp_metadata.json").read_text())).set_index("column")
    check("metadata agrees with exported columns", set(metadata.index) == set(raw) - {"Unnamed: 0"}, True)
    check("MIC columns", int(metadata.output_kind.eq("MIC_uM").sum()), 51)
    check("finite positive MIC exports", bool(np.isfinite(raw[metadata.index[metadata.output_kind.eq("MIC_uM")]].stack()).all() and raw[metadata.index[metadata.output_kind.eq("MIC_uM")]].stack().gt(0).all()), True)
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence", verify_integrity=True)
    check("exact source coverage", set(source.index) == set(raw.index), True)
    # Independently reconstruct source labels using pinned local upstream tables.
    # Upstream checkout's append file is a Git LFS pointer. Verify and use the
    # already materialized local cache rather than parsing that pointer as data.
    append_cache = ROOT / "scratchpad/marlys_source.tsv"
    append_pointer = (ROOT / "scratchpad/DAMP/data/data_damp/df_marlys_append.tsv").read_text()
    append_oid = next(line.split("sha256:")[1] for line in append_pointer.splitlines() if line.startswith("oid sha256:"))
    check("appended source cache matches upstream LFS hash", hashlib.sha256(append_cache.read_bytes()).hexdigest(), append_oid)
    upstream = pd.concat([
        pd.read_csv(ROOT / "scratchpad/DAMP/data/df_damp_190973x4.tsv", sep="\t"),
        pd.read_csv(append_cache, sep="\t"),
    ], ignore_index=True).set_index("sequence", verify_integrity=True)
    check("source assignments match pinned tables", bool(source.source.eq(upstream.source.reindex(source.index)).all()), True)
    seq = raw.index.to_series()
    check("canonical amino acid alphabet", bool(seq.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+").all()), True)
    f = pd.DataFrame({"source": source.source.reindex(raw.index), "length": seq.str.len()}, index=raw.index)
    f["charge_density"] = (seq.str.count("[KR]") - seq.str.count("[DE]")) / f.length
    f["cys_fraction"] = seq.str.count("C") / f.length
    f["hydrophobic_fraction"] = seq.str.count("[AILMFWVY]") / f.length
    f["entropy_bits"] = [-sum((c / len(s)) * np.log2(c / len(s)) for c in Counter(s).values()) for s in seq]
    values = raw[list(COLS.values())].rename(columns={v: k for k, v in COLS.items()}).dropna()
    z = np.log2(values)
    f = f.loc[z.index].copy()
    primary = z[PRIMARY]
    spread = primary.max(axis=1) - primary.min(axis=1)
    low, high = spread.quantile([.05, .95])
    f["group"] = np.select([spread.le(low), spread.ge(high)], ["low", "high"], default="middle")
    f["log2_spread"] = spread
    f["median_log2_prediction"] = primary.median(axis=1)
    rank = primary.rank(pct=True)
    f["rank_spread"] = rank.max(axis=1) - rank.min(axis=1)
    centered = primary - primary.median()
    centered_spread = centered.max(axis=1) - centered.min(axis=1)
    threshold = pd.read_csv(ROOT / "analysis_outputs/high_low_thresholds.csv").set_index("population").loc["all"]
    check("complete cases", len(z), threshold.n)
    check("low tail cutoff log2", low, threshold.low_max_log2)
    check("high tail cutoff log2", high, threshold.high_min_log2)
    old_profiles = pd.read_csv(ROOT / "analysis_outputs/high_low_profiles.csv").query("population == 'all'").set_index("group")
    for group, block in f.groupby("group"):
        old = old_profiles.loc[group]
        check(f"{group} n", len(block), old.n)
        for feature in ["length", "charge_density", "cys_fraction", "hydrophobic_fraction", "entropy_bits"]:
            check(f"{group} {feature} mean", block[feature].mean(), old[feature + "_mean"])
    tolerance = pd.read_csv(ROOT / "analysis_outputs/mic_variability_tolerance.csv")
    for _, row in tolerance.iterrows():
        check(f"raw fraction within {row.maximum_fold_spread:g}fold", spread.le(np.log2(row.maximum_fold_spread)).mean(), row.raw_fraction)
        check(f"centered fraction within {row.maximum_fold_spread:g}fold", centered_spread.le(np.log2(row.maximum_fold_spread)).mean(), row.centered_fraction)
    old_high = pd.read_csv(ROOT / "analysis_outputs/high_low_sequences_all.csv").set_index("sequence")
    check("all saved group memberships", bool(f.group.eq(old_high.group.reindex(f.index)).all()), True)
    for old_name, new_name in OLD_NAMES.items():
        if old_name in old_high:
            check(f"saved values {new_name}", bool(np.allclose(values[new_name], old_high[old_name].reindex(values.index))), True)
    # Check archived matching definitions, uniqueness, and mean feature differences.
    matches = pd.read_csv(ROOT / "analysis_outputs/all_high_vs_rest_matches.csv")
    a = f.loc[matches.high_sequence].reset_index(drop=True)
    b = f.loc[matches.control_sequence].reset_index(drop=True)
    check("matching source and exact length", bool((a.source.eq(b.source) & a.length.eq(b.length)).all()), True)
    check("matching predicted MIC bin", bool(np.floor(a.median_log2_prediction).eq(np.floor(b.median_log2_prediction)).all()), True)
    check("matching no replacement", matches.control_sequence.is_unique and matches.high_sequence.is_unique, True)
    matched = pd.read_csv(ROOT / "analysis_outputs/high_low_matched.csv").query("population == 'all' and contrast == 'high_vs_rest'").set_index("feature")
    for feature in matched.index:
        check(f"matched {feature}", (a[feature] - b[feature]).mean(), matched.at[feature, "mean_paired_difference"])
    # Quantify the finite-length constraint on composition entropy without treating
    # residues as iid draws or substituting a correction with unverified assumptions.
    f["entropy_source_length_residual"] = f.entropy_bits - f.groupby(["source", "length"]).entropy_bits.transform("mean")
    entropy_rows = []
    for group, block in f.groupby("group"):
        entropy_rows.append({"group": group, "n": len(block), "mean_entropy_bits": block.entropy_bits.mean(),
                             "mean_source_exact_length_residual_bits": block.entropy_source_length_residual.mean()})
    save(pd.DataFrame(entropy_rows), "entropy_length_sensitivity")

    # Terminal-extension components across the full input, not just extreme rows.
    parent = {s: s for s in raw.index}
    def root(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s
    pairs = []
    for s in raw.index:
        if len(s) < 9:
            continue
        for t in {s[1:], s[:-1]}:
            if t in parent:
                parent[root(s)] = root(t)
                pairs.append((s, t))
    check("terminal extension pairs", len(pairs), 92233)
    check("terminal participant sequences", len({s for pair in pairs for s in pair}), 87340)
    families = pd.Series({s: root(s) for s in raw.index}, name="terminal_component")
    f["terminal_component"] = families.reindex(f.index)
    # One fixed representative per component, selected without looking at predictions.
    hashes = pd.Series({s: hashlib.sha256(s.encode()).hexdigest() for s in z.index})
    representatives = hashes.groupby(f.terminal_component).idxmin().to_numpy()
    restriction = f.length.between(8, 39) & f.cys_fraction.eq(0)
    populations = {"all": z, "cysteine_free_8_to_39": z.loc[restriction],
                   "one_terminal_component_representative": z.loc[representatives]}
    pair_rows = []
    shortlist_rows = []
    for population, zz in populations.items():
        pair_rows.extend(pair_metrics(zz, population))
        zz = zz.sort_index()  # deterministic tie breaks
        for fraction in [.005, .01, .05, .1]:
            k = int(np.ceil(len(zz) * fraction))
            selected = {model: zz[model].nsmallest(k).index for model in zz}
            for left, right in combinations(zz, 2):
                overlap = len(selected[left].intersection(selected[right]))
                record = {"population": population, "fraction": fraction, "left": left, "right": right,
                          "n": len(zz), "k": k, "overlap_count": overlap, "overlap_fraction": overlap / k,
                          "jaccard": overlap / (2 * k - overlap), "expected_unconditional_count": k * k / len(zz)}
                # Uniform independent selection within strata, preserving each model's
                # selected counts there. This is a reference, not a significance test.
                for label, grouping in [("source", ["source"]), ("source_length", ["source", "length"])]:
                    universe = f.loc[zz.index].groupby(grouping).size()
                    nl = f.loc[selected[left]].groupby(grouping).size().reindex(universe.index, fill_value=0)
                    nr = f.loc[selected[right]].groupby(grouping).size().reindex(universe.index, fill_value=0)
                    expected = float((nl * nr / universe).sum())
                    record[f"expected_{label}_count"] = expected
                    record[f"observed_over_{label}_reference"] = overlap / expected if expected else np.nan
                shortlist_rows.append(record)
                if population == "all":
                    old = pd.read_csv(ROOT / "analysis_outputs/shortlist_sensitivity.csv")
                    old = old[(old.fraction == fraction) & old.left.map(LONG_OLD).eq(left) & old.right.map(LONG_OLD).eq(right)]
                    check(f"shortlist {fraction} {left} / {right}", overlap, int(old.iloc[0].overlap_count))
    pair_table = pd.DataFrame(pair_rows)
    shortlist = pd.DataFrame(shortlist_rows)
    save(pair_table, "scale_and_ranking")
    save(shortlist, "shortlists_with_selection_references")

    # Exact decomposition of pooled log prediction covariance into between-source
    # covariance of means plus mean within-source covariance, using population ddof=0.
    covariance_rows = []
    for left, right in combinations(z, 2):
        aa, bb = z[left], z[right]
        total = float(np.mean((aa - aa.mean()) * (bb - bb.mean())))
        amean = aa.groupby(f.source).transform("mean")
        bmean = bb.groupby(f.source).transform("mean")
        between = float(np.mean((amean - aa.mean()) * (bmean - bb.mean())))
        within = float(np.mean((aa - amean) * (bb - bmean)))
        check(f"covariance identity {left}/{right}", within + between, total)
        covariance_rows.append({"left": left, "right": right, "total_log2_covariance": total,
                                "between_source_covariance": between, "within_source_covariance": within,
                                "between_fraction_of_covariance": between / total})
    save(pd.DataFrame(covariance_rows), "source_covariance_decomposition")

    # Agreement tails: interpret each model's relative position as well as its value.
    profiles = []
    for group in ["low", "middle", "high"]:
        idx = f.index[f.group.eq(group)]
        for model in PRIMARY:
            profiles.append({"group": group, "model": model, "n": len(idx),
                             "median_prediction_uM": values.loc[idx, model].median(),
                             "median_within_model_percentile": rank.loc[idx, model].median(),
                             "group_median_rank_spread": f.loc[idx, "rank_spread"].median(),
                             "group_median_centered_fold_spread": np.exp2(centered_spread.loc[idx].median())})
    save(pd.DataFrame(profiles), "agreement_group_interpretation")
    # Cross-definition tail transition matrix, rows normalize to raw-tail group sizes.
    definitions = {"raw_log2_spread": spread, "median_centered_spread": centered_spread,
                   "within_model_percentile_spread": f.rank_spread}
    transition_rows = []
    feature_rows = []
    for name, score in definitions.items():
        lo, hi = score.quantile([.05, .95])
        groups = pd.Series(np.select([score.le(lo), score.ge(hi)], ["low", "high"], default="middle"), index=f.index)
        for raw_group in ["low", "middle", "high"]:
            subset = groups.loc[f.group.eq(raw_group)]
            for alternative in ["low", "middle", "high"]:
                transition_rows.append({"definition": name, "raw_group": raw_group, "alternative_group": alternative,
                                        "n": int(subset.eq(alternative).sum()), "fraction": subset.eq(alternative).mean()})
        for group in ["low", "high"]:
            b = f.loc[groups.eq(group)]
            feature_rows.append({"definition": name, "group": group, "n": len(b),
                                 "median_length": b.length.median(), "cysteine_present_fraction": b.cys_fraction.gt(0).mean(),
                                 "mean_charge_density": b.charge_density.mean(), "mean_entropy_bits": b.entropy_bits.mean()})
    save(pd.DataFrame(transition_rows), "agreement_definition_transitions")
    save(pd.DataFrame(feature_rows), "sequence_profiles_by_agreement_definition")

    # A deliberately stringent local contrast: both models must change at least
    # twofold in opposite directions for an extension to count as a strong reversal.
    pair_frame = pd.DataFrame(pairs, columns=["longer", "shorter"])
    pair_frame = pair_frame[pair_frame.longer.isin(z.index) & pair_frame.shorter.isin(z.index)].copy()
    deltas = z.loc[pair_frame.longer].to_numpy() - z.loc[pair_frame.shorter].to_numpy()
    delta = pd.DataFrame(deltas, columns=z.columns)
    terminal_rows = []
    for left, right in combinations(z, 2):
        da, db = delta[left], delta[right]
        strong_both = da.abs().ge(1) & db.abs().ge(1)
        opposite = da.mul(db).lt(0)
        terminal_rows.append({"left": left, "right": right, "n_pairs": len(delta),
                              "n_both_change_twofold": int(strong_both.sum()),
                              "n_opposite_twofold_changes": int((strong_both & opposite).sum()),
                              "fraction_all_pairs_opposite_twofold": (strong_both & opposite).mean(),
                              "fraction_of_strong_pairs_opposite": opposite.loc[strong_both].mean(),
                              "median_extension_effect_left_log2": da.median(),
                              "median_extension_effect_right_log2": db.median()})
    save(pd.DataFrame(terminal_rows), "terminal_extension_prediction_changes")
    example_mask = delta["APEX 1.1"].ge(1) & delta["AMPredictor"].le(-1)
    examples = pair_frame.iloc[np.flatnonzero(example_mask)].copy()
    for model in z:
        examples[model + "_longer_uM"] = values.loc[examples.longer, model].to_numpy()
        examples[model + "_shorter_uM"] = values.loc[examples.shorter, model].to_numpy()
    examples["source_longer"] = f.loc[examples.longer, "source"].to_numpy()
    examples["source_shorter"] = f.loc[examples.shorter, "source"].to_numpy()
    save(examples, "opposed_extension_examples")

    # Falsification-oriented check: extreme disagreement is not automatically an
    # APEX inactive-label ceiling phenomenon. 500 is a proximity diagnostic only.
    ceiling_rows = []
    for group, b in f.groupby("group"):
        ceiling_rows.append({"group": group, "n": len(b),
                             "fraction_APEX_1_1_ge_500": values.loc[b.index, "APEX 1.1"].ge(500).mean(),
                             "fraction_APEX_1_1_ge_512": values.loc[b.index, "APEX 1.1"].ge(512).mean()})
    save(pd.DataFrame(ceiling_rows), "inactive_label_proximity")

    # Source-dependent offsets: learn a centering constant on all OTHER sources,
    # apply to held-out source. This measures transfer of a prediction discrepancy,
    # not predictive calibration or generalization against experimental labels.
    transfer = []
    for name, idx in f.groupby("source").groups.items():
        hold = primary.loc[idx]
        train = primary.drop(index=idx)
        test_centered = hold - train.median()
        local = hold - hold.median()
        transfer.append({"source": name, "n": len(hold),
                         "raw_median_fold_spread": np.exp2((hold.max(axis=1) - hold.min(axis=1)).median()),
                         "other_source_centered_median_fold_spread": np.exp2((test_centered.max(axis=1) - test_centered.min(axis=1)).median()),
                         "within_source_centered_median_fold_spread": np.exp2((local.max(axis=1) - local.min(axis=1)).median())})
    save(pd.DataFrame(transfer), "source_offset_transfer")

    # Independently recompute the small experimental diagnostic from the curated
    # labels and check exposure. This does not endorse independence of these labels.
    e = pd.read_csv(ROOT / "analysis_outputs/experimental_ecoli_matched.csv").set_index("sequence")
    metrics = pd.read_csv(ROOT / "analysis_outputs/experimental_metrics.csv").set_index("model")
    db = pd.read_csv(ROOT / "scratchpad/replay/APEXDB.csv").set_index("Sequence", verify_integrity=True)
    check("numeric experimental sample", len(e), 13)
    check("APEXDB exact exposure of numeric matches", bool(e.index.isin(db.index).all()), True)
    for old, new in OLD_NAMES.items():
        error = np.log2(raw.loc[e.index, COLS[new]]) - np.log2(e.MIC_uM)
        check(f"diagnostic MAE {new}", error.abs().mean(), metrics.at[old, "MAE_log2"])
    replay = pd.read_csv(ROOT / "analysis_outputs/deepamp_replay.csv")
    token_max = float(replay.loc[replay.padding.eq("token_padding"), "log2_ratio"].abs().max())
    check("DeepAMP replay actual maximum < 6.2e-6", token_max < 6.2e-6, True)
    check("DeepAMP saved replay contains 160 token-padding values", int(replay.padding.eq("token_padding").sum()), 160)
    save(pd.DataFrame(checks), "independent_checks")

    # Figures with concrete explanations, shared denominators and exportable PDF.
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "pdf.fonttype": 42})
    colors = {"APEX 1.1": "#225b78", "MBC-Attention": "#c16936", "AMPredictor": "#287a65"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.9), layout="constrained")
    grid = np.linspace(0, 10, 501)
    for score, label, color in [(spread, "Exported MIC predictions", "#b04a41"),
                                 (centered_spread, "After removing each model's median", "#225b78")]:
        axes[0].plot(2 ** grid, np.searchsorted(np.sort(score), grid, side="right") / len(score) * 100, label=label, color=color, lw=2.3)
    axes[0].set(xscale="log", xticks=[2,4,8,16,32,128,512], xlim=(1,1024), ylim=(0,100),
                xlabel="Largest / smallest prediction (fold)", ylabel="Sequences within this spread (%)",
                title="A. Scale offsets account for much of the spread")
    axes[0].set_xticklabels([2,4,8,16,32,128,512])
    axes[0].legend(fontsize=8, loc="lower right")
    selected = shortlist[(shortlist.population == "all") & shortlist.fraction.eq(.01) & shortlist.left.isin(PRIMARY)]
    labels = [r.left + "\nvs " + r.right for r in selected.itertuples()]
    axes[1].barh(labels, selected.overlap_fraction.to_numpy() * 100, color=["#225b78", "#287a65", "#c16936"])
    for i, r in enumerate(selected.itertuples()):
        axes[1].text(r.overlap_fraction * 100 + .13, i, f"{r.overlap_count:,} / {r.k:,} ({r.overlap_fraction:.1%})", va="center", fontsize=9)
    axes[1].axvline(1, color="#777777", ls="--", lw=1)
    axes[1].set(xlim=(0,8.3), xlabel="Lowest-MIC 1% shortlist shared (%)",
                title="B. Removing offsets cannot change the shortlist")
    axes[1].invert_yaxis()
    fig.suptitle("The models differ in both numerical scale and candidate priority", fontsize=14, weight="bold")
    figsave(fig, "figure_1_scale_and_priority")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    profile = pd.DataFrame(profiles)
    positions = np.arange(3)
    for i, model in enumerate(PRIMARY):
        yy = profile[profile.model.eq(model)].set_index("group").loc[["low", "middle", "high"], "median_within_model_percentile"] * 100
        axes[0].bar(positions + (i-1)*.24, yy, width=.23, label=model, color=colors[model])
    axes[0].set(xticks=positions, xticklabels=["Low MIC spread\n(bottom 5%)", "Middle 90%", "High MIC spread\n(top 5%)"],
                ylim=(0,100), ylabel="Median within-model percentile (lower = more potent)",
                title="A. Similar MICs can have very different relative ranks")
    axes[0].legend(fontsize=8)
    for group, color in [("low", "#225b78"), ("middle", "#a5a9ad"), ("high", "#b04a41")]:
        score = f.loc[f.group.eq(group), "rank_spread"].sort_values()
        axes[1].plot(score.to_numpy()*100, np.arange(1,len(score)+1)/len(score)*100, color=color, lw=2, label=group)
    axes[1].set(xlim=(0,100), ylim=(0,100), xlabel="Within-sequence percentile spread (percentage points)",
                ylabel="Sequences within this rank spread (%)", title="B. Both raw MIC tails show large rank discrepancies")
    axes[1].legend(title="Raw MIC-spread group", fontsize=8)
    fig.suptitle("Numerical agreement should not be treated as screening confidence", fontsize=14, weight="bold")
    figsave(fig, "figure_2_what_agreement_means")

    sources = sorted(f.source.unique())
    correlation = []
    for name in sources:
        block = z.loc[f.source.eq(name)]
        row = {"source": name, "n": len(block)}
        for left, right in combinations(PRIMARY, 2):
            row[left + " / " + right] = block[left].corr(block[right], method="spearman")
        correlation.append(row)
    source_corr = pd.DataFrame(correlation).set_index("source")
    save(source_corr.reset_index(), "within_source_correlations")
    fig, axes = plt.subplots(1,2,figsize=(12,6),layout="constrained", gridspec_kw={"width_ratios":[1.2,1]})
    matrix = source_corr.drop(columns="n")
    im = axes[0].imshow(matrix, vmin=-1, vmax=1, cmap="RdBu", aspect="auto")
    axes[0].set(yticks=np.arange(len(matrix)), yticklabels=[s.replace("_", " ")+f" (n={int(source_corr.at[s,'n']):,})" for s in matrix.index],
                xticks=np.arange(3), xticklabels=["APEX /\nMBC", "APEX /\nAMPredictor", "MBC /\nAMPredictor"],
                title="A. Ranking agreement depends on the collection")
    for i in range(len(matrix)):
        for j in range(3):
            axes[0].text(j,i,f"{matrix.iloc[i,j]:.2f}",ha="center",va="center",color="white" if abs(matrix.iloc[i,j])>.6 else "#17262d")
    fig.colorbar(im, ax=axes[0], shrink=.65, label="Spearman correlation")
    for population, color, label in [("all", "#225b78", "All complete cases"),
                                      ("cysteine_free_8_to_39", "#c16936", "Cysteine-free; 8–39 residues"),
                                      ("one_terminal_component_representative", "#287a65", "One per terminal component")]:
        rows = shortlist[(shortlist.population == population) & shortlist.fraction.eq(.01) & shortlist.left.isin(PRIMARY)]
        offset = {"all":-.24,"cysteine_free_8_to_39":0,"one_terminal_component_representative":.24}[population]
        axes[1].barh(np.arange(3)+offset,rows.overlap_fraction*100,height=.23,color=color,label=label)
    axes[1].set(yticks=np.arange(3), yticklabels=["APEX / MBC", "APEX / AMPredictor", "MBC / AMPredictor"],
                xlabel="Lowest-MIC 1% shortlist shared (%)", title="B. Small overlaps persist in sensitivity analyses")
    axes[1].invert_yaxis();axes[1].legend(fontsize=8,loc="upper center",bbox_to_anchor=(.5,-.12))
    fig.suptitle("Collection composition and sequence redundancy matter",fontsize=14,weight="bold")
    figsave(fig,"figure_3_sources_and_robustness")

    summary = {"input_sha256": checksum, "python": platform.python_version(), "numpy": np.__version__,
               "pandas": pd.__version__, "scipy": scipy.__version__, "n_sequences": len(raw),
               "n_complete": len(z), "n_restricted": int(restriction.sum()),
               "n_terminal_representatives": len(representatives),
               "terminal_components_spanning_sources": int(f.groupby("terminal_component").source.nunique().gt(1).sum()),
               "independent_checks": len(checks), "checks_passed": sum(c["passed"] for c in checks),
               "check_discrepancies": [c for c in checks if not c["passed"]],
               "median_raw_fold_spread": float(np.exp2(spread.median())),
               "median_centered_fold_spread": float(np.exp2(centered_spread.median())),
               "original_low_group_median_percentile_spread": float(f.loc[f.group.eq("low"), "rank_spread"].median()),
               "original_high_group_median_percentile_spread": float(f.loc[f.group.eq("high"), "rank_spread"].median()),
               "DeepAMP_saved_replay_max_absolute_log2_difference": token_max,
               "reporting_correction": "Original <6e-6 bound is slightly too small; corrected bound is <6.2e-6 log2 units.",
               "new_inference_performed_by_this_script": False,
               "inferential_scope": "Exploratory finite-collection prediction behavior; no experimental accuracy or population significance."}
    (OUT / "review_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("\nPrimary 1% shortlists:\n",shortlist[(shortlist.population == "all") & shortlist.fraction.eq(.01) & shortlist.left.isin(PRIMARY)].to_string(index=False))
    print("\nAgreement group profiles:\n",profile.to_string(index=False))
    print("\nTerminal extension contrasts:\n",pd.DataFrame(terminal_rows).to_string(index=False))


if __name__ == "__main__":
    main()
