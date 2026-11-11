"""Compare all eight DAMP sources along the final report's main dimensions.

Run: python analysis_review/source_stratified_analysis.py
Uses saved predictions and observed sequence links; no new model inference.
New results are isolated in outputs/source_comparisons.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from classifier_sequence_analysis import CLASS, matched_contrasts
from plot_shortlist_upset import MIC

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/source_comparisons"
ORDER = ["Venomics", "marlys_8to50_novel", "Archaeomics", "Human_Proteomics", "Extinctomics", "DRAMP_APD_DBAASP", "NEURIPS_AMP_CHALLENGE", "AMPDIFF_synth50k"]
DEEP = {target: [f"deep_amp_{model}_{suffix}_unspecified_species_and_strain" for model in ["cnn", "lstm"]]
        for target, suffix in [("E. coli", "neg_Gram-negative"), ("B. subtilis", "pos_Gram-positive")]}


def save(data, name):
    pd.DataFrame(data).to_csv(OUT / (name + ".csv"), index=False)


def shortlists(predictions, frame):
    rows, ranks = [], []
    for source, f in frame.groupby("source"):
        p = predictions.loc[f.index]
        populations = {"every_sequence": p,
                       "C_free_length8_39": p.loc[f.contains_C.eq(0) & f.length.between(8, 39)],
                       "group_medians": p.groupby(f.component).median()}
        for seed in range(20):
            hashes = pd.Series([hashlib.sha256(f"{seed}|{s}".encode()).hexdigest() for s in p.index], index=p.index)
            chosen = f.loc[hashes.sort_values().index].drop_duplicates("component").index
            populations[f"representative_{seed:02d}"] = p.loc[chosen]
        for population, data in populations.items():
            for left, right in combinations(MIC, 2):
                if population in ["every_sequence", "group_medians", "C_free_length8_39"]:
                    ranks.append({"source": source, "population": population, "n": len(data), "left": left, "right": right,
                                  "spearman": data[left].corr(data[right], method="spearman")})
            for fraction in ([.01, .02, .05, .1] if population == "every_sequence" else [.05]):
                k = int(np.ceil(len(data) * fraction))
                selected = {m: set(data[m].nsmallest(k).index) for m in MIC}
                triple = len(set.intersection(*selected.values()))
                for left, right in combinations(MIC, 2):
                    shared = len(selected[left] & selected[right])
                    rows.append({"source": source, "population": population, "fraction": fraction, "n_entries": len(data),
                                 "k_each": k, "left": left, "right": right, "shared_n": shared,
                                 "shared_fraction": shared / k, "triple_shared_n": triple,
                                 "boundary_tie_left": int(data[left].eq(data[left].nsmallest(k).iloc[-1]).sum()) > 1,
                                 "boundary_tie_right": int(data[right].eq(data[right].nsmallest(k).iloc[-1]).sum()) > 1})
    table = pd.DataFrame(rows)
    save(table, "shortlist_comparisons")
    save(ranks, "rank_correlations")
    reps = table[table.population.str.startswith("representative_")]
    summary = reps.groupby(["source", "left", "right"]).agg(n_entries=("n_entries", "first"), k_each=("k_each", "first"),
                runs=("shared_fraction", "size"), minimum=("shared_fraction", "min"), median=("shared_fraction", "median"), maximum=("shared_fraction", "max")).reset_index()
    save(summary, "representative_shortlist_summary")
    return table


def matched_decisions(raw, source, component):
    specs = {**DEEP, "APEX E. coli versions": ["apex_1_0_E._coli_ATCC11775", "apex_1_1_E._coli_ATCC_11775"]}
    rows = []
    for target, cols in specs.items():
        paired = raw[cols].dropna()
        for label, ix in source.reindex(paired.index).groupby(source.reindex(paired.index)).groups.items():
            p = paired.loc[ix]
            mapped = component.reindex(ix).notna()
            group_ids = component.reindex(ix[mapped])
            weights = 1 / group_ids.map(group_ids.value_counts())
            for cutoff in [8, 16, 32]:
                bits = p.le(cutoff)
                both = bits.all(axis=1); either = bits.any(axis=1); disagree = bits.iloc[:, 0] != bits.iloc[:, 1]
                strong = ((p.iloc[:, 0].le(cutoff / 2) & p.iloc[:, 1].gt(cutoff * 2)) |
                          (p.iloc[:, 1].le(cutoff / 2) & p.iloc[:, 0].gt(cutoff * 2)))
                weighted_xor = float(np.sum(weights * disagree.reindex(weights.index)))
                weighted_any = float(np.sum(weights * either.reindex(weights.index)))
                rows.append({"source": label, "target": target, "cutoff_uM": cutoff, "n_paired": len(p),
                             "both_pass_n": int(both.sum()), "both_fail_n": int((~either).sum()), "either_pass_n": int(either.sum()),
                             "disagree_n": int(disagree.sum()), "overall_disagreement": disagree.mean(),
                             "selected_pool_disagreement": disagree.sum() / either.sum() if either.any() else np.nan,
                             "strong_margin_n": int(strong.sum()), "strong_margin_fraction": strong.mean(),
                             "n_component_mapped": int(mapped.sum()), "n_groups": group_ids.nunique(),
                             "group_weighted_overall_disagreement": weighted_xor / weights.sum(),
                             "group_weighted_selected_disagreement": weighted_xor / weighted_any if weighted_any else np.nan})
    save(rows, "matched_target_decisions")
    return pd.DataFrame(rows)


def order_responses(raw, pairs):
    rows = []
    for source, p in pairs.groupby("source"):
        a = raw.loc[p.sequence_left]; b = raw.loc[p.sequence_right]
        for target, columns in DEEP.items():
            va = a[columns].to_numpy(); vb = b[columns].to_numpy()
            delta = np.abs(np.log2(vb / va))
            populations = {"all_observed_pairs": np.ones(len(p), bool),
                           "no_group_reuse": p.greedy_component_representative.to_numpy(dtype=bool),
                           "all_predictions_le128": np.all((va <= 128) & (vb <= 128), axis=1)}
            for population, keep in populations.items():
                for j, model in enumerate(["CNN", "LSTM"]):
                    rows.append({"source": source, "target": target, "population": population, "model": model,
                                 "n_pairs": int(keep.sum()), "twofold_n": int((delta[keep, j] >= 1).sum()),
                                 "twofold_fraction": np.mean(delta[keep, j] >= 1) if keep.any() else np.nan})
    save(rows, "sequence_order_responses")
    result = pd.DataFrame(rows)
    balanced = result.groupby(["target", "population", "model"]).apply(
        lambda g: pd.Series({"n_sources": g.twofold_fraction.notna().sum(),
                             "n_pairs": g.n_pairs.sum(),
                             "equal_source_fraction": g.twofold_fraction.mean(),
                             "pair_weighted_fraction": g.twofold_n.sum() / g.n_pairs.sum()})).reset_index()
    save(balanced, "source_balanced_order")
    return result


def classifier_context(frame):
    rows = []
    for source, f in frame.groupby("source"):
        for population, weights in [("every_sequence", pd.Series(1., index=f.index)),
                                    ("equal_group", 1 / f.groupby("component").component.transform("size"))]:
            result = matched_contrasts(f, weights, population)
            rows.extend({"source": source, **row} for row in result if row["left"] == "All negative" and row["right"] == "All positive" and row["feature"] in ["KR_fraction", "DE_fraction"])
    save(rows, "classifier_length_matched_profiles")
    # The catalog's observed fractions remain the same as the previous analysis;
    # reference probabilities are now computed within the primary MIC cohort.
    pairs = pd.read_csv(ROOT / "analysis_review/outputs/classifier_sequences/one_substitution_pairs.csv.gz", usecols=["sequence_from", "sequence_to", "source", "component"])
    counts = frame.groupby(["source", "length", "group"]).group.transform("size")
    total = frame.groupby(["source", "length"]).group.transform("size")
    probability = (counts - 1) / (total - 1).replace(0, np.nan)
    rows = []
    for source, p in pairs.groupby("source"):
        left = p.sequence_from; right = p.sequence_to
        same = frame.group.reindex(left).to_numpy() == frame.group.reindex(right).to_numpy()
        reference = (probability.reindex(left).to_numpy() + probability.reindex(right).to_numpy()) / 2
        for population, weight in [("every_pair", np.ones(len(p))), ("equal_pair_group", 1 / p.component.map(p.component.value_counts()).to_numpy())]:
            rows.append({"source": source, "population": population, "n_pairs": len(p), "n_groups": p.component.nunique(),
                         "same_group_fraction": np.average(same, weights=weight), "source_length_reference": np.average(reference, weights=weight)})
    save(rows, "neighbor_agreement")


def figures(shortlists, decisions, order):
    pairs = list(combinations(MIC, 2)); pair_labels = ["APEX / MBC", "APEX / AMPredictor", "MBC / AMPredictor"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), layout="constrained")
    for ax, pop, title in zip(axes, ["every_sequence", "group_medians"], ["A. Every sequence", "B. Each related group's median prediction"]):
        d = shortlists[shortlists.population.eq(pop) & shortlists.fraction.eq(.05)]
        matrix = np.array([[d[d.source.eq(s) & d.left.eq(a) & d.right.eq(b)].shared_fraction.iloc[0] * 100 for a, b in pairs] for s in ORDER])
        im = ax.imshow(matrix, vmin=0, vmax=60, cmap="Blues", aspect="auto")
        for i in range(8):
            for j in range(3): ax.text(j, i, f"{matrix[i,j]:.1f}%", ha="center", va="center", color="white" if matrix[i,j] > 38 else "#20313c")
        sizes = d.groupby("source").first()
        labels = [f"{s} (n={int(sizes.loc[s,'n_entries']):,}; k={int(sizes.loc[s,'k_each']):,})" for s in ORDER]
        ax.set(yticks=range(8), yticklabels=labels, xticks=range(3), xticklabels=pair_labels, title=title)
        ax.tick_params(axis="x", labelrotation=15)
    fig.colorbar(im, ax=axes, label="Shared count / either top-5% list size", shrink=.8)
    fig.suptitle("Source-specific shortlists: do the same peptides or groups get selected?", weight="bold", fontsize=14)
    figsave(fig, "source_shortlist_overlap")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.0), layout="constrained")
    d = decisions[decisions.target.eq("E. coli") & decisions.cutoff_uM.eq(16)].set_index("source").loc[ORDER]
    for ax, field, title, nfield in [(axes[0], "overall_disagreement", "A. All paired predictions", "n_paired"), (axes[1], "selected_pool_disagreement", "B. Peptides passing CNN or LSTM", "either_pass_n")]:
        ax.barh(range(8), d[field] * 100, color=["#225b78" if s == "Venomics" else "#a5b8c3" for s in ORDER])
        for i, val in enumerate(d[field] * 100): ax.text(val + 1, i, f"{val:.1f}%", va="center", fontsize=9)
        ax.set(yticks=range(8), yticklabels=[f"{s} (n={int(d.loc[s,nfield]):,})" for s in ORDER], xlim=(0, 80), xlabel="% with exactly one model passing", title=title)
        ax.invert_yaxis()
    fig.suptitle("DeepAMP E. coli decisions at 16 µM: the denominator changes the impression", weight="bold", fontsize=14)
    figsave(fig, "source_matched_decision_denominators")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6), layout="constrained")
    for ax, target in zip(axes, DEEP):
        data = order[order.target.eq(target)]
        for j, (model, color) in enumerate([("CNN", "#225b78"), ("LSTM", "#bd6a34")]):
            d = data[data.model.eq(model) & data.population.eq("all_observed_pairs")].set_index("source").loc[ORDER]
            ax.barh(np.arange(8) + (j - .5) * .32, d.twofold_fraction * 100, height=.30, color=color, alpha=.75, label=model + ": all pairs")
            r = data[data.model.eq(model) & data.population.eq("no_group_reuse")].set_index("source").loc[ORDER]
            ax.scatter(r.twofold_fraction * 100, np.arange(8) + (j - .5) * .32, color=color, marker="D", s=30, edgecolors="black", linewidth=.5, label=model + ": no group reuse")
        n = data.groupby(["source", "population"]).n_pairs.first()
        labels = [f"{s} ({n.loc[s,'all_observed_pairs']}/{n.loc[s,'no_group_reuse']} pairs)" for s in ORDER]
        ax.set(yticks=range(8), yticklabels=labels, xlim=(0, 100), xlabel="% with at least a twofold order response", title=target)
        ax.invert_yaxis(); ax.legend(loc="upper center", bbox_to_anchor=(.5, -.12), ncol=2, fontsize=8)
    fig.suptitle("Sequence order: source patterns and pair coverage differ", weight="bold", fontsize=14)
    figsave(fig, "source_order_responses")


def figsave(fig, name):
    for ext in ["png", "pdf", "svg"]: fig.savefig(OUT / (name + "." + ext), dpi=190, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    paths = [ROOT / "df_damp_195625x57.tsv", ROOT / "scratchpad/source_sequences.tsv", ROOT / "analysis_review/outputs/classifier_sequences/sequence_groups.csv.gz", ROOT / "analysis_review/outputs/sanity_checks/broader_sequence_components.csv", ROOT / "analysis_review/outputs/screening_decision_trends/DeepAMP_composition_controlled_pairs.csv"]
    raw = pd.read_csv(paths[0], sep="\t").set_index("sequence", verify_integrity=True)
    source = pd.read_csv(paths[1], sep="\t").set_index("sequence", verify_integrity=True).source
    component = pd.read_csv(paths[3]).set_index("sequence", verify_integrity=True).component
    predictions = raw[list(MIC.values())].rename(columns={c: m for m, c in MIC.items()}).dropna().sort_index()
    frame = pd.read_csv(paths[2]).set_index("sequence", verify_integrity=True).loc[predictions.index]
    assert len(frame) == 194947 and frame.source.notna().all() and frame.component.notna().all()
    coverage = []
    for label in ORDER:
        ix = source.index[source.eq(label)]
        coverage.append({"source": label, "full_n": len(ix), "primary_MIC_n": int(frame.source.eq(label).sum()),
                         "classifier_complete_n": int(raw.loc[ix, list(CLASS.values())].notna().all(axis=1).sum()),
                         "DeepAMP_paired_n": int(raw.loc[ix, DEEP["E. coli"]].notna().all(axis=1).sum())})
    save(coverage, "coverage")
    shortlist = shortlists(predictions, frame)
    decisions = matched_decisions(raw, source, component)
    classifier_context(frame)
    pairs = pd.read_csv(paths[4])
    order = order_responses(raw, pairs)
    # Cached primary source summaries are explicitly verified, then linked by the report.
    cached = pd.read_csv(ROOT / "analysis_review/outputs/venomics/within_source_shortlists.csv")
    join = shortlist[shortlist.population.eq("every_sequence")].merge(cached, on=["source", "fraction", "left", "right"], suffixes=("_new", "_old"))
    old_decisions = pd.read_csv(ROOT / "analysis_review/outputs/screening_decision_trends/matched_target_by_source.csv")
    old_deep = old_decisions[old_decisions.comparison.eq("DeepAMP_same_Gram_group")].copy()
    old_deep["target"] = old_deep.target.map({"Gram-negative": "E. coli", "Gram-positive": "B. subtilis"})
    joined = decisions[decisions.cutoff_uM.eq(16) & decisions.target.isin(DEEP)].merge(old_deep, on=["source", "target"], suffixes=("_new", "_old"))
    checks = {"all_eight_sources_primary_194947": int(sum(r["primary_MIC_n"] for r in coverage)) == 194947,
              "all_full_rows_195625": sum(r["full_n"] for r in coverage) == 195625,
              "96_prior_shortlist_intersections_reproduced": len(join) == 96 and bool(join.shared_n_new.eq(join.shared_n_old).all()) and bool(np.allclose(join.shared_fraction, join.fraction_of_either_list)),
              "16_DeepAMP_source_decision_counts_reproduced": len(joined) == 16 and bool(joined.disagree_n.eq(joined.decision_disagreement_n).all()),
              "matched_cells_partition_each_cohort": bool((decisions.both_pass_n + decisions.both_fail_n + decisions.disagree_n).eq(decisions.n_paired).all()),
              "order_source_counts_partition_5306": int(order[(order.target.eq("E. coli")) & order.model.eq("CNN") & order.population.eq("all_observed_pairs")].n_pairs.sum()) == 5306,
              "no_group_reuse_counts_partition_2197": int(order[(order.target.eq("E. coli")) & order.model.eq("CNN") & order.population.eq("no_group_reuse")].n_pairs.sum()) == 2197,
              "all_order_pair_ingredients_identical": all(sorted(a) == sorted(b) for a, b in zip(pairs.sequence_left, pairs.sequence_right)),
              "all_order_pairs_same_source": bool(source.reindex(pairs.sequence_left).to_numpy().tolist() == source.reindex(pairs.sequence_right).to_numpy().tolist()),
              "all_primary_source_indices_match": bool(frame.source.eq(source.reindex(frame.index)).all())}
    assert all(checks.values()), checks
    save([{"check": k, "passed": bool(v)} for k, v in checks.items()], "checks")
    manifest = {"scope": "Exploratory descriptive source stratification of ranking, concentration decisions, matched targets, classifier sequence profiles, neighbor decisions and composition-controlled order responses. No source or model accuracy ranking.",
                "source_order": ORDER, "primary_n": 194947, "matched_n": "All paired outputs, separately reported per source and target; mapped-only subsets for group weights.",
                "representatives": "20 deterministic SHA256(seed|sequence) choices per source-component; ranges describe choices, not uncertainty intervals",
                "neighbor_reference": "Same source and exact length in the primary MIC cohort, excluding the anchor; pair weighting and equal pair-component weighting both reported.",
                "software": {"pandas": pd.__version__, "numpy": np.__version__},
                "input_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                "checks": {k: bool(v) for k, v in checks.items()}, "all_checks_passed": all(checks.values())}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    figures(shortlist, decisions, order)
    print(f"Source comparison: all {len(checks)} checks passed; outputs in {OUT.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
