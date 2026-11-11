"""Exact shortlist and 16 µM cutoff UpSets, plus matched-sequence rank correlations.

Run: python analysis_review/plot_shortlist_upset.py
All budgets use the same 194,947 complete cases and ceil(fraction * N) entries.
Exclusive UpSet intersections differ from inclusive pairwise overlaps.
Spearman describes ordering, not concentration agreement or biological accuracy.
"""
from pathlib import Path
from itertools import product, combinations
import argparse
import hashlib
import json
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter, MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/shortlist_overlap"
MIC = {"APEX 1.1": "apex_1_1_E._coli_ATCC_11775",
       "MBC-Attention": "mbc_attention_E._coli_unspecified_strain",
       "AMPredictor": "ampredictor_E._coli_unspecified_strain"}
FRACTIONS = [.01, .02, .05]


def save_figure(fig, directory, name):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(directory / (name + "." + ext), dpi=200, bbox_inches="tight")
    plt.close(fig)


def draw_upset(counts, order, k, n, fraction, directory, *, selection_sizes=None,
               title=None, subtitle=None, stem="shortlist_upset"):
    fig = plt.figure(figsize=(10.8, 6.3), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[2.8, 7.6], height_ratios=[3.6, 1.6], hspace=.03, wspace=.08)
    bars = fig.add_subplot(grid[0, 1]); matrix = fig.add_subplot(grid[1, 1], sharex=bars)
    sizes = fig.add_subplot(grid[1, 0], sharey=matrix)
    xx = np.arange(len(order)); heights = [counts[p] for p in order]
    colors = ["#225b78" if sum(p) == 1 else ("#bd6a34" if sum(p) == 2 else "#287a65") for p in order]
    bars.bar(xx, heights, width=.67, color=colors)
    for j, count in enumerate(heights):
        bars.text(j, count + max(heights) * .02, f"{count:,}", ha="center", va="bottom", weight="bold", fontsize=11)
    bars.set(ylabel="Peptides in this exact intersection", ylim=(0, max(heights) * 1.15))
    bars.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    bars.tick_params(axis="x", bottom=False, labelbottom=False)
    for j, pattern in enumerate(order):
        yy = np.flatnonzero(pattern)
        matrix.scatter([j] * 3, range(3), s=46, color="#dce3e7", zorder=1)
        if len(yy) > 1: matrix.plot([j, j], [min(yy), max(yy)], color=colors[j], lw=2.3, zorder=2)
        matrix.scatter([j] * len(yy), yy, s=65, color=colors[j], zorder=3)
    matrix.set(ylim=(2.6, -.6), xlim=(-.6, len(order)-.4), yticks=range(3), xticks=[])
    matrix.tick_params(axis="y", left=False, labelleft=False)
    matrix.set_xlabel("Filled dots identify the models selecting that group of peptides", labelpad=12, fontsize=10)
    for spine in matrix.spines.values(): spine.set_visible(False)
    selection_sizes = [k] * 3 if selection_sizes is None else selection_sizes
    sizes.barh(range(3), selection_sizes, color="#a5b8c3", height=.46)
    sizes.set(yticks=range(3), yticklabels=list(MIC), xlim=(max(selection_sizes) * 1.13, 0), xlabel="Selected peptides" if k is None else "Shortlist size")
    sizes.xaxis.set_major_locator(MaxNLocator(3, integer=True))
    sizes.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
    sizes.tick_params(axis="y", length=0)
    for y, size in enumerate(selection_sizes): sizes.text(size * .955, y, f"{size:,}", va="center", ha="left", fontsize=9)
    sizes.spines["left"].set_visible(False)
    fig.suptitle(title or f"Overlap of the lowest-MIC top {fraction:.0%} lists", fontsize=17, weight="bold")
    bars.set_title(subtitle or f"{k:,} peptides per model • same {n:,}-sequence collection • exclusive intersections", fontsize=10, pad=16)
    save_figure(fig, directory, stem)


def cutoff_overlap(v, cutoff=16):
    """Exact cutoff memberships, including joint failures in the audit table."""
    directory = OUT / "cutoff_16uM"
    directory.mkdir(parents=True, exist_ok=True)
    membership = v.le(cutoff)
    membership.to_csv(directory / "cutoff_membership.csv.gz", compression="gzip")
    patterns = list(product([0, 1], repeat=3))
    counts = {p: int(np.all(membership.to_numpy() == np.array(p), axis=1).sum()) for p in patterns}
    order = sorted([p for p in patterns if any(p)], key=lambda p: (sum(p), -counts[p], p))
    pd.DataFrame([{"pattern": "".join(map(str, p)), "exclusive_count": counts[p],
                   "selected_by_n_models": sum(p), **{m: int(p[j]) for j, m in enumerate(MIC)}}
                  for p in [(0, 0, 0)] + order]).to_csv(directory / "exclusive_intersections.csv", index=False)
    selected = {m: set(v.index[v[m] <= cutoff]) for m in MIC}
    universe = set(v.index)
    sizes = [len(selected[m]) for m in MIC]
    union = set.union(*selected.values())
    checks = {
        "common_cohort_194947": len(v) == 194947,
        "all_eight_cells_partition_full_cohort": sum(counts.values()) == len(v),
        "seven_displayed_cells_partition_any_pass": sum(counts[p] for p in order) == len(union),
        "cell_counts_reproduce_each_model_total": all(sum(counts[p] for p in patterns if p[j]) == sizes[j] for j in range(3)),
        "independent_set_operations_reproduce_all_eight_cells": all(
            len(set.intersection(*(selected[m] if p[j] else universe - selected[m] for j, m in enumerate(MIC)))) == counts[p]
            for p in patterns),
        "totals_match_reported_cutoff_analysis": sizes == [1362, 100676, 126192] and len(union) == 144361 and counts[(1, 1, 1)] == 897,
        "no_missing_or_nonpositive_predictions": bool(np.isfinite(v).all().all() and v.gt(0).all().all()),
    }
    record = {"cutoff_uM": cutoff, "rule": "exported predicted MIC <= cutoff",
              "input_sha256": hashlib.sha256((ROOT / "df_damp_195625x57.tsv").read_bytes()).hexdigest(),
              "n_common": len(v), "selected_counts": dict(zip(MIC, sizes)), "any_pass": len(union),
              "all_pass": counts[(1, 1, 1)], "none_pass": counts[(0, 0, 0)],
              "checks": checks, "all_checks_passed": all(checks.values()),
              "scope": "Same approximately matched E. coli outputs as the shortlist comparison. No clinical breakpoint or accuracy inference."}
    (directory / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    assert record["all_checks_passed"], record
    draw_upset(counts, order, None, len(v), None, directory, selection_sizes=sizes,
               title="Overlap of peptides with predicted MIC ≤16 µM",
               subtitle=f"Same {len(v):,} sequences • {len(union):,} pass any model • {counts[(0, 0, 0)]:,} pass none",
               stem="cutoff_upset")
    return record


def shortlist(v, fraction, directory):
    directory.mkdir(parents=True, exist_ok=True)
    k = int(np.ceil(fraction * len(v)))
    selected = {m: set(v[m].nsmallest(k).index) for m in MIC}
    union = sorted(set.union(*selected.values()))
    membership = pd.DataFrame({m: [s in selected[m] for s in union] for m in MIC}, index=pd.Index(union, name="sequence"))
    membership.to_csv(directory / "shortlist_membership.csv")
    patterns = [p for p in product([0, 1], repeat=3) if any(p)]
    counts = {p: int(np.all(membership.to_numpy() == np.array(p), axis=1).sum()) for p in patterns}
    order = sorted(patterns, key=lambda p: (sum(p), -counts[p], p))
    rows = [{"pattern": "".join(map(str, p)), "exclusive_count": counts[p], "selected_by_n_models": sum(p),
             **{m: int(p[j]) for j, m in enumerate(MIC)}} for p in order]
    pd.DataFrame(rows).to_csv(directory / "exclusive_intersections.csv", index=False)
    pairwise = [{"left": left, "right": right, "k": k,
                 "intersection_including_triple": len(selected[left] & selected[right]),
                 "fraction_of_either_shortlist": len(selected[left] & selected[right]) / k}
                for left, right in combinations(MIC, 2)]
    pd.DataFrame(pairwise).to_csv(directory / "inclusive_pairwise_overlaps.csv", index=False)
    checks = {
        "common_cohort_194947": len(v) == 194947,
        ("each_shortlist_has_1950" if fraction == .01 else "each_shortlist_has_declared_size"): all(len(s) == k for s in selected.values()),
        "exclusive_bars_partition_union": sum(counts.values()) == len(union),
        "exclusive_memberships_reconstruct_total_selections": sum(sum(p) * counts[p] for p in counts) == 3 * k,
        "exclusive_memberships_reconstruct_each_shortlist": all(sum(counts[p] for p in counts if p[j]) == k for j in range(3)),
        "independent_set_algebra_reproduces_all_exclusive_cells": all(
            len(set.intersection(*(selected[m] if p[j] else set(union) - selected[m] for j, m in enumerate(MIC)))) == counts[p]
            for p in counts),
        "no_tied_selection_boundary": all(v[m].eq(v[m].nsmallest(k).iloc[-1]).sum() == 1 for m in MIC),
    }
    if fraction == .01: checks["pairwise_counts_reproduce_previous_report"] = [r["intersection_including_triple"] for r in pairwise] == [42, 33, 107]
    manifest = {"input_sha256": hashlib.sha256((ROOT / "df_damp_195625x57.tsv").read_bytes()).hexdigest(),
                "n_common": len(v), "shortlist_size_each": k, "union_n": len(union), "all_three_n": counts[(1, 1, 1)],
                "one_model_only_n": sum(counts[p] for p in counts if sum(p) == 1),
                "shortlist_cutoffs_uM": {m: float(v[m].nsmallest(k).iloc[-1]) for m in MIC},
                "checks": checks, "all_checks_passed": all(checks.values()),
                "scope": "Exact exclusive overlaps of model-selected candidates; E. coli endpoints are approximately matched. No accuracy or confidence inference."}
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    assert all(checks.values()), manifest
    draw_upset(counts, order, k, len(v), fraction, directory)
    summary = {"top_percent": int(fraction * 100), "n_common": len(v), "k": k, "union_n": len(union),
               "one_model_only_n": manifest["one_model_only_n"],
               "two_models_only_n": sum(counts[p] for p in counts if sum(p) == 2), "all_three_n": counts[(1, 1, 1)]}
    for row in pairwise:
        row["top_percent"] = int(fraction * 100)
        row["jaccard_intersection_over_union"] = row["intersection_including_triple"] / (2 * k - row["intersection_including_triple"])
    return summary, pairwise, checks


def correlations(v):
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence", verify_integrity=True).source.reindex(v.index)
    assert source.notna().all()
    populations = {"All shared sequences": v.index, **source.groupby(source).groups}
    rows = []
    rank_checks = []
    for population, index in populations.items():
        data = v.loc[index]; corr = data.corr(method="spearman")
        independent = data.rank(method="average").corr(method="pearson")
        rank_checks.append(bool(np.allclose(corr, independent, atol=1e-12, rtol=0)))
        for left, right in combinations(MIC, 2):
            rows.append({"population": population, "left": left, "right": right, "n": len(data), "spearman": corr.loc[left, right]})
    table = pd.DataFrame(rows); table.to_csv(OUT / "spearman_correlations.csv", index=False)
    cached = pd.read_csv(ROOT / "analysis_review/outputs/scale_and_ranking.csv").query("population == 'all'")
    global_rows = table[table.population.eq("All shared sequences")]
    assert all(np.isclose(row.spearman, cached[(cached.left.eq(row.left)) & cached.right.eq(row.right)].spearman.iloc[0], atol=1e-12)
               for row in global_rows.itertuples())
    table["pair"] = table.left + " / " + table.right
    matrix = table.pivot(index="population", columns="pair", values="spearman").reindex(index=populations, columns=global_rows.left + " / " + global_rows.right)
    fig, ax = plt.subplots(figsize=(10.5, 6.6), layout="constrained")
    plot = ax.imshow(matrix, cmap="RdBu", vmin=-1, vmax=1, aspect="auto")
    for i in range(len(matrix)):
        for j in range(3): ax.text(j, i, f"{matrix.iloc[i, j]:.3f}", ha="center", va="center", fontsize=11)
    n = table.groupby("population").n.first()
    labels = [f"{name.replace('_', ' ')} (n = {n[name]:,})" for name in matrix.index]
    ax.set(xticks=range(3), xticklabels=["APEX / MBC", "APEX / AMPredictor", "MBC / AMPredictor"],
           yticks=range(len(matrix)), yticklabels=labels, title="Rank agreement across models, overall and within each source")
    ax.get_yticklabels()[0].set_weight("bold")
    fig.colorbar(plot, ax=ax, label="Spearman rank correlation (−1 to +1)")
    save_figure(fig, OUT, "spearman_correlations")
    return {"average_tied_rank_Pearson_matches_Spearman_in_every_population": all(rank_checks),
            "global_correlations_match_previous_independent_review": True,
            "all_eight_sources_and_full_cohort_reported": len(populations) == 9}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff-only", action="store_true", help="Build only the 16 µM overlap, preserving prior result files")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    path = ROOT / "df_damp_195625x57.tsv"
    v = pd.read_csv(path, sep="\t", usecols=["sequence"] + list(MIC.values())).set_index("sequence", verify_integrity=True)
    v = v.rename(columns={col: name for name, col in MIC.items()}).dropna().sort_index()
    assert np.isfinite(v).all().all() and v.gt(0).all().all()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    cutoff_record = cutoff_overlap(v)
    if args.cutoff_only:
        print(json.dumps(cutoff_record, indent=2))
        return
    summaries = []; pairs = []; budget_checks = {}
    for fraction in FRACTIONS:
        directory = OUT if fraction == .01 else OUT / f"top_{int(fraction * 100)}pct"
        summary, pairwise, checks = shortlist(v, fraction, directory)
        summaries.append(summary); pairs.extend(pairwise); budget_checks[str(fraction)] = checks
    pd.DataFrame(summaries).to_csv(OUT / "budget_summary.csv", index=False)
    pd.DataFrame(pairs).to_csv(OUT / "budget_pairwise_overlaps.csv", index=False)
    rank_checks = correlations(v)
    record = {"input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "n_common": len(v),
              "fractions": FRACTIONS, "rounding": "ceil(fraction * N)", "budget_checks": budget_checks,
              "correlation_checks": rank_checks, "all_checks_passed": all(rank_checks.values()) and all(all(c.values()) for c in budget_checks.values()),
              "scope": "Spearman uses the same sequences within each stated population, not ranks restricted to a selected shortlist. No biological accuracy or significance claim."}
    (OUT / "budget_extension_manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    assert record["all_checks_passed"], record
    if __package__:
        from .report_figures import build_decision_panels
    else:
        from report_figures import build_decision_panels
    build_decision_panels()
    print(json.dumps({"budgets": summaries, "all_checks_passed": record["all_checks_passed"]}, indent=2))


if __name__ == "__main__":
    main()
