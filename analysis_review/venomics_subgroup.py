"""Venomics versus the other DAMP sources, using saved predictions only.

Run: python analysis_review/venomics_subgroup.py
New outputs are isolated under outputs/venomics. No fitting or model inference.
"""
from pathlib import Path
from itertools import combinations, product
import hashlib
import json
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from classifier_sequence_analysis import CLASS
from plot_shortlist_upset import MIC, draw_upset

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/venomics"
FEATURES = ["KR_fraction", "DE_fraction", "contains_C", "hydrophobic_fraction"]
OUTCOMES = ["all_positive", "mixed", "all_negative", "APEX_pass16", "MBC_pass16", "AMP_pass16"]


def save(rows, name):
    pd.DataFrame(rows).to_csv(OUT / (name + ".csv"), index=False)


def compare(z, keys, label, weights=None):
    """Standardize other sequences to Venomics' distribution in shared cells.

    Each cell is weighted by its Venomics mass, rather than minimum cell count.
    No causal interpretation: the retained population and support are explicit.
    """
    z = z.copy()
    z["population"] = np.where(z.source.eq("Venomics"), "Venomics", "Other sources")
    z["mass"] = 1. if weights is None else weights.reindex(z.index)
    grouped = z.groupby(keys + ["population"], observed=True)
    mass = grouped.mass.sum().unstack("population").reindex(columns=["Venomics", "Other sources"]).dropna()
    rows = []
    for feature in OUTCOMES + FEATURES:
        weighted = z[feature] * z.mass
        means = weighted.groupby([z[k] for k in keys] + [z.population], observed=True).sum().unstack("population").reindex(index=mass.index, columns=mass.columns) / mass
        ven = float(np.average(means.Venomics, weights=mass.Venomics))
        other = float(np.average(means["Other sources"], weights=mass.Venomics))
        rows.append({"comparison": label, "metric": feature, "Venomics_fraction": ven, "other_fraction": other,
                     "difference_percentage_points": 100 * (ven - other), "common_cells": len(mass),
                     "Venomics_mass_in_common_cells": mass.Venomics.sum(),
                     "Venomics_total_mass": z.loc[z.source.eq("Venomics"), "mass"].sum()})
    return rows


def source_summary(frame, predictions):
    rows, correlations, cutoff = [], [], []
    for source, z in frame.groupby("source"):
        sizes = z.component.value_counts()
        weights = 1 / z.component.map(sizes)
        row = {"source": source, "n": len(z), "n_groups": len(sizes), "median_length": z.length.median(),
               "largest_group_n": sizes.max(), "top10_groups_fraction": sizes.nlargest(10).sum() / len(z),
               "training_composition_proxy_n": int((z.contains_C.eq(0) & z.length.between(8, 39)).sum())}
        for feature in FEATURES + OUTCOMES:
            row[feature] = z[feature].mean()
            row[feature + "_equal_group"] = np.average(z[feature], weights=weights)
        rows.append(row)
        p = predictions.loc[z.index]
        for left, right in combinations(MIC, 2):
            correlations.append({"source": source, "left": left, "right": right, "n": len(z),
                                 "spearman": p[left].corr(p[right], method="spearman")})
        for model in MIC:
            for threshold in [4, 8, 16, 32, 64]:
                cutoff.append({"source": source, "model": model, "cutoff_uM": threshold, "n": len(z),
                               "pass_n": int(p[model].le(threshold).sum()), "pass_fraction": p[model].le(threshold).mean()})
    save(rows, "source_context")
    save(correlations, "source_correlations")
    save(cutoff, "source_cutoffs")
    return pd.DataFrame(rows)


def overlaps(predictions):
    patterns = list(product([0, 1], repeat=3))
    rows, cells = [], []
    for source, p in predictions.groupby("source"):
        p = p.drop(columns="source").sort_index()
        for fraction in [.01, .02, .05, .10]:
            k = int(np.ceil(fraction * len(p)))
            bits = pd.DataFrame({model: p.index.isin(p[model].nsmallest(k).index) for model in MIC}, index=p.index)
            for left, right in combinations(MIC, 2):
                shared = int((bits[left] & bits[right]).sum())
                rows.append({"source": source, "fraction": fraction, "k": k, "n": len(p), "left": left, "right": right,
                             "shared_n": shared, "fraction_of_either_list": shared / k,
                             "boundary_tied_left": int(p[left].eq(p[left].nsmallest(k).iloc[-1]).sum()) > 1,
                             "boundary_tied_right": int(p[right].eq(p[right].nsmallest(k).iloc[-1]).sum()) > 1})
            if source == "Venomics":
                counts = {b: int(np.all(bits.to_numpy() == np.array(b), axis=1).sum()) for b in patterns}
                cells.extend({"rule": f"top_{fraction:.0%}", "pattern": "".join(map(str, b)), "n": counts[b]} for b in patterns)
                assert all(sum(counts[b] for b in patterns if b[j]) == k for j in range(3))
                if fraction == .05:
                    selected = [b for b in patterns if any(b)]
                    order = sorted(selected, key=lambda b: (sum(b), -counts[b], b))
                    draw_upset(counts, order, k, len(p), fraction, OUT,
                               title="Within Venomics: overlap of the top 5% MIC lists", stem="venomics_top5_upset")
                    bits.to_csv(OUT / "venomics_top5_membership.csv", index_label="sequence")
                if fraction == .01:
                    assert sum(counts.values()) == len(p)
        if source == "Venomics":
            bits = p.le(16)
            counts = {b: int(np.all(bits.to_numpy() == np.array(b), axis=1).sum()) for b in patterns}
            cells.extend({"rule": "MIC_le_16", "pattern": "".join(map(str, b)), "n": counts[b]} for b in patterns)
            bits.to_csv(OUT / "venomics_cutoff_membership.csv", index_label="sequence")
    save(rows, "within_source_shortlists")
    save(cells, "venomics_exact_overlaps")


def shortlist_sensitivities(frame, predictions):
    """Check the Venomics ranking result without counting every variant equally."""
    v = frame[frame.source.eq("Venomics")]
    p = predictions.loc[v.index]
    populations = {"all_Venomics": p, "no_C_length_8_39": p.loc[v.contains_C.eq(0) & v.length.between(8, 39)],
                   "group_medians": p.groupby(v.component).median()}
    for seed in range(20):
        rank = pd.Series([hashlib.sha256(f"{seed}|{s}".encode()).hexdigest() for s in v.index], index=v.index)
        chosen = v.loc[rank.sort_values().index].drop_duplicates("component").index
        populations[f"one_per_group_seed_{seed:02d}"] = p.loc[chosen]
    rows = []
    for population, data in populations.items():
        k = int(np.ceil(.05 * len(data)))
        selected = {m: set(data[m].nsmallest(k).index) for m in MIC}
        triple = len(set.intersection(*selected.values()))
        for left, right in combinations(MIC, 2):
            shared = len(selected[left] & selected[right])
            rows.append({"population": population, "n_entries": len(data), "k_each": k, "left": left, "right": right,
                         "shared_n": shared, "shared_fraction": shared / k, "triple_shared_n": triple})
    save(rows, "shortlist_sensitivities")


def make_plots(context, comparisons):
    order = context.sort_values("all_positive", ascending=False).source.tolist()
    d = context.set_index("source").loc[order]
    labels = [f"{s.replace('_', ' ')} (n={int(d.loc[s, 'n']):,})" for s in order]
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2), layout="constrained")
    bottom = np.zeros(len(d))
    for feature, label, color in [("all_negative", "All negative", "#8a9da8"), ("mixed", "Mixed", "#bd6a34"), ("all_positive", "All positive", "#287a65")]:
        values = d[feature].to_numpy() * 100
        axes[0].barh(range(len(d)), values, left=bottom, label=label, color=color)
        if feature == "all_positive":
            for j, value in enumerate(values): axes[0].text(101, j, f"{value:.1f}%", va="center", fontsize=9)
        bottom += values
    axes[0].set(yticks=range(len(d)), yticklabels=labels, xlim=(0, 118), xlabel="% of that source; right labels = all positive", title="A. Classifier groups (positive score ≥0.5)")
    axes[0].invert_yaxis()
    axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.12), ncol=3, fontsize=9)
    values = d[["APEX_pass16", "MBC_pass16", "AMP_pass16"]].to_numpy() * 100
    im = axes[1].imshow(values, vmin=0, vmax=100, cmap="Blues", aspect="auto")
    for i in range(len(d)):
        for j in range(3): axes[1].text(j, i, f"{values[i,j]:.2f}%", ha="center", va="center", color="white" if values[i,j] > 65 else "#20313c")
    axes[1].set(yticks=range(len(d)), yticklabels=order, xticks=range(3), xticklabels=["APEX 1.1", "MBC", "AMPredictor"], title="B. Predicted MIC ≤16 µM")
    fig.colorbar(im, ax=axes[1], label="% passing within that source", shrink=.8)
    fig.suptitle("Venomics: frequent classifier consensus, different MIC screening choices", weight="bold", fontsize=15)
    figsave(fig, "venomics_screening_context")
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 5.6), layout="constrained")
    for ax, feature, title in zip(axes, ["KR_fraction", "contains_C", "top10_groups_fraction"], ["A. Mean lysine/arginine (K/R) content", "B. Peptides containing cysteine (C)", "C. Sequences in the ten largest groups"]):
        ax.barh(range(len(d)), d[feature] * 100, color=["#225b78" if s == "Venomics" else "#a5b8c3" for s in order])
        if feature != "top10_groups_fraction":
            ax.scatter(d[feature + "_equal_group"] * 100, range(len(d)), color="#bd6a34", s=35, label="Equal related-group weight", zorder=3)
            ax.legend(loc="upper center", bbox_to_anchor=(.5, -.12), fontsize=8)
        for i, value in enumerate(d[feature] * 100): ax.text(value + 1, i, f"{value:.1f}", va="center", fontsize=8)
        ax.set(yticks=range(len(d)), yticklabels=order if ax is axes[0] else [], xlabel="%", title=title, xlim=(0, max(100 if feature == "contains_C" else 65, d[feature].max() * 115)))
        ax.invert_yaxis()
    fig.suptitle("Sequence ingredients and sampling differ between DAMP sources", weight="bold", fontsize=15)
    figsave(fig, "venomics_sequence_context")
    fig, ax = plt.subplots(figsize=(10.5, 5), layout="constrained")
    metrics = ["all_positive", "APEX_pass16", "MBC_pass16", "AMP_pass16"]
    colors = ["#225b78", "#bd6a34", "#287a65", "#8b548a", "#707070"]
    labels = {"raw": "Every sequence", "equal_group": "Every related group", "length_C": "Same length + cysteine presence", "length_C_KR": "Also similar K/R content", "equal_group_length_C_KR": "Group weights + sequence matching"}
    plotted = comparisons[comparisons.comparison.isin(labels)]
    for j, (name, block) in enumerate(plotted.groupby("comparison", sort=False)):
        data = block.set_index("metric").loc[metrics]
        ax.scatter(data.difference_percentage_points, np.arange(4) + (j - 2) * .13, label=labels[name], color=colors[j], s=45)
    ax.axvline(0, color="#8a9da8", lw=.9)
    ax.set(yticks=range(4), yticklabels=["All four classifiers positive", "APEX 1.1: MIC ≤16", "MBC: MIC ≤16", "AMPredictor: MIC ≤16"], xlabel="Venomics − other sources (percentage points)", title="Which differences survive changes in weighting and sequence context?")
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2, fontsize=9)
    figsave(fig, "venomics_comparison_controls")


def figsave(fig, name):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / (name + "." + ext), dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [ROOT / "df_damp_195625x57.tsv", ROOT / "scratchpad/source_sequences.tsv", ROOT / "analysis_review/outputs/classifier_sequences/sequence_groups.csv.gz"]
    raw = pd.read_csv(paths[0], sep="\t").set_index("sequence", verify_integrity=True)
    source_table = pd.read_csv(paths[1], sep="\t").set_index("sequence", verify_integrity=True)
    predictions = raw[list(MIC.values())].rename(columns={c: m for m, c in MIC.items()}).dropna().sort_index()
    frame = pd.read_csv(paths[2]).set_index("sequence", verify_integrity=True).loc[predictions.index].copy()
    scores = raw.loc[frame.index, list(CLASS.values())]
    assert len(frame) == 194947 and frame.source.eq("Venomics").sum() == 4618
    assert frame.source.eq(source_table.source.reindex(frame.index)).all() and frame.component.notna().all()
    assert scores.notna().all().all() and np.isfinite(predictions).all().all() and predictions.gt(0).all().all()
    for name, count in [("all_positive", 4), ("all_negative", 0)]: frame[name] = scores.ge(.5).sum(axis=1).eq(count).astype(float)
    frame["mixed"] = 1 - frame.all_positive - frame.all_negative
    for name, model in [("APEX_pass16", "APEX 1.1"), ("MBC_pass16", "MBC-Attention"), ("AMP_pass16", "AMPredictor")]: frame[name] = predictions[model].le(16).astype(float)
    frame["KR_bin"] = np.floor(frame.KR_fraction * 10 + 1e-10).astype(int)
    context = source_summary(frame, predictions)
    weights = 1 / frame.groupby(["source", "component"]).component.transform("size")
    comparisons = []
    for label, w in [("raw", pd.Series(1., index=frame.index)), ("equal_group", weights)]:
        for feature in OUTCOMES + FEATURES:
            v = frame.source.eq("Venomics")
            a = np.average(frame.loc[v, feature], weights=w[v]); b = np.average(frame.loc[~v, feature], weights=w[~v])
            comparisons.append({"comparison": label, "metric": feature, "Venomics_fraction": a, "other_fraction": b,
                                "difference_percentage_points": 100 * (a - b), "common_cells": np.nan,
                                "Venomics_mass_in_common_cells": w[v].sum(), "Venomics_total_mass": w[v].sum()})
    comparisons += compare(frame, ["length", "contains_C"], "length_C")
    comparisons += compare(frame, ["length", "contains_C", "KR_bin"], "length_C_KR")
    comparisons += compare(frame, ["length", "contains_C", "KR_bin"], "equal_group_length_C_KR", weights=weights)
    bin_sensitivity = []
    for width in [.05, .20]:
        alternative = frame.copy()
        alternative["KR_bin"] = np.floor(alternative.KR_fraction / width + 1e-10).astype(int)
        bin_sensitivity.extend({"KR_bin_width": width, **r} for r in compare(alternative, ["length", "contains_C", "KR_bin"], f"equal_group_KR_width_{width}", weights=weights))
    save(bin_sensitivity, "matching_bin_sensitivity")
    source_controls = []
    for source in sorted(set(frame.source) - {"Venomics"}):
        selected = frame[frame.source.isin(["Venomics", source])]
        source_controls.extend({"other_source": source, **r} for r in compare(selected, ["length", "contains_C", "KR_bin"], "equal_group_length_C_KR", weights=weights))
    save(source_controls, "one_source_at_a_time_controls")
    comparisons = pd.DataFrame(comparisons)
    save(comparisons, "controlled_comparisons")
    # Equal-source comparator keeps each other collection equally influential.
    equal_sources = frame.groupby("source")[OUTCOMES + FEATURES].mean()
    save(equal_sources.reset_index(), "source_mean_metrics")
    save([{"metric": f, "Venomics_fraction": equal_sources.loc["Venomics", f],
           "other_source_mean_fraction": equal_sources.drop(index="Venomics")[f].mean()} for f in OUTCOMES + FEATURES], "equal_source_reference")
    p = predictions.assign(source=frame.source)
    overlaps(p)
    shortlist_sensitivities(frame, predictions)
    # Compare both classifiers and MIC cutoffs under alternate score rules.
    threshold_rows = []
    for threshold in [.3, .5, .7]:
        votes = scores.ge(threshold).sum(axis=1)
        for source, ix in frame.groupby("source").groups.items():
            v = votes.loc[ix]
            threshold_rows.append({"source": source, "threshold": threshold, "n": len(v),
                                   "all_positive_n": int(v.eq(4).sum()), "mixed_n": int(v.between(1, 3).sum()), "all_negative_n": int(v.eq(0).sum())})
    save(threshold_rows, "classifier_rule_sensitivity")
    # Selected-source context: eleven original APEX targets used by the mining paper.
    old_targets = ["E._coli_ATCC11775", "P._aeruginosa_PAO1", "P._aeruginosa_PA14", "S._aureus_ATCC12600", "E._coli_AIG221", "E._coli_AIG222", "K._pneumoniae_ATCC13883", "A._baumannii_ATCC19606", "S._aureus_(ATCC_BAA-1556)_-_MRSA", "vancomycin-resistant_E._faecalis_ATCC700802", "vancomycin-resistant_E._faecium_ATCC700221"]
    old_median = raw.loc[frame.index, ["apex_1_0_" + t for t in old_targets]].median(axis=1)
    save([{"population": label, "n": int(mask.sum()), "median_of_panel_medians_uM": old_median[mask].median(),
           "panel_median_le32_n": int(old_median[mask].le(32).sum()), "panel_median_le32_01_n": int(old_median[mask].le(32.01).sum()),
           "maximum_panel_median_uM": old_median[mask].max()} for label, mask in [("Venomics", frame.source.eq("Venomics")), ("Other sources", frame.source.ne("Venomics"))]], "original_selection_context")
    frame.loc[frame.source.eq("Venomics"), ["source", "length", "component", *FEATURES, *OUTCOMES]].join(predictions).to_csv(OUT / "venomics_sequence_annotations.csv", index_label="sequence")
    checks = {"primary_cohort_194947": len(frame) == 194947, "Venomics_all_4618_retained": frame.source.eq("Venomics").sum() == source_table.source.eq("Venomics").sum() == 4618,
              "classifier_groups_partition_each_source": bool(np.allclose(context[["all_positive", "mixed", "all_negative"]].sum(axis=1), 1)),
              "Venomics_401_groups": int(context.set_index("source").loc["Venomics", "n_groups"]) == 401,
              "Venomics_cutoffs_reproduce_32_3987_4212": [int(frame.loc[frame.source.eq("Venomics"), f].sum()) for f in OUTCOMES[3:]] == [32, 3987, 4212],
              "matched_length_C_features_exact": all(abs(row.difference_percentage_points) < 1e-10 for row in comparisons[(comparisons.comparison.eq("length_C")) & comparisons.metric.eq("contains_C")].itertuples()),
              "all_controlled_fractions_bounded": bool(comparisons[["Venomics_fraction", "other_fraction"]].ge(0).all().all() and comparisons[["Venomics_fraction", "other_fraction"]].le(1).all().all())}
    assert all(checks.values()), checks
    save([{"check": k, "passed": bool(v)} for k, v in checks.items()], "checks")
    manifest = {"scope": "Exploratory source-specific model behavior, not venom-wide prevalence or accuracy; Venomics is already prediction-selected.",
                "n_common": len(frame), "n_Venomics": 4618, "n_other": len(frame)-4618, "models": MIC,
                "classifier_rule": "score >= 0.5; 0.3 and 0.7 sensitivity; AMP-CLIP binary",
                "matched_weights": "Venomics sequence counts within common exact-length/C-presence strata; extra KR bins are floor(KR_fraction*10), 10 percentage points wide. The comparator uses pooled non-Venomics sequences within each cell. Controls are descriptive, not causal.",
                "equal_group_weights": "1/(observed component size within source); each represented source-group has equal total influence, not each source",
                "source_origin": {"paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12254355/", "import_code": "scratchpad/1_DATA_PROBE.py:61-98", "source_table": "scratchpad/source_sequences.tsv"},
                "input_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}, "checks": {k: bool(v) for k,v in checks.items()}, "all_checks_passed": all(checks.values())}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    make_plots(context, comparisons)
    print(context.to_string(index=False))
    print(comparisons[comparisons.metric.isin(OUTCOMES)].to_string(index=False))


if __name__ == "__main__":
    main()
