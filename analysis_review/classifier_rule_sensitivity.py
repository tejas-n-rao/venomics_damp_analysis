"""Check classifier groups at permissive (0.3), baseline (0.5), strict (0.7) rules.

Run: python analysis_review/classifier_rule_sensitivity.py
Uses the same complete-score cohort and observed neighbors; no new inference.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from classifier_sequence_analysis import CLASS, GROUPS, group_labels, matched_contrasts

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis_review/outputs/classifier_sequences"
OUT = BASE / "rule_sensitivity"
RULES = [.3, .5, .7]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw_path = ROOT / "df_damp_195625x57.tsv"
    frame_path = BASE / "sequence_groups.csv.gz"
    pair_path = BASE / "one_substitution_pairs.csv.gz"
    scores = pd.read_csv(raw_path, sep="\t", usecols=["sequence", *CLASS.values()]).set_index("sequence").dropna().sort_index()
    frame = pd.read_csv(frame_path).set_index("sequence").reindex(scores.index)
    pairs = pd.read_csv(pair_path, usecols=["sequence_from", "sequence_to"])
    assert len(scores) == 195620 and len(pairs) == 41253
    assert frame.source.notna().all() and scores.ge(0).all().all() and scores.le(1).all().all()
    assert set(scores[CLASS["AMP-CLIP"]].unique()) <= {0, 1}
    # Independently check the two features used in the plotted matched contrasts.
    seq = scores.index.to_series()
    assert np.allclose(frame.KR_fraction, seq.str.count("[KR]") / seq.str.len())
    assert np.allclose(frame.DE_fraction, seq.str.count("[DE]") / seq.str.len())
    votes_rows, contrasts, neighbors, transitions = [], [], [], []
    checks = {}
    cached = pd.read_csv(BASE / "score_rule_sensitivity.csv")
    previous_votes = None
    for threshold in RULES:
        bits = scores.ge(threshold)
        votes = bits.sum(axis=1)
        z = frame.copy()
        z["group"] = group_labels(votes)
        mixed_n = int(votes.between(1, 3).sum())
        for vote in range(5):
            count = int(votes.eq(vote).sum())
            votes_rows.append({"threshold": threshold, "positive_votes": vote, "n": count,
                               "fraction_of_all": count / len(z),
                               "fraction_of_mixed": count / mixed_n if 0 < vote < 4 else np.nan})
        contrasts.extend({"threshold": threshold, **r} for r in matched_contrasts(z, pd.Series(1., index=z.index), "all_complete"))
        for group in GROUPS:
            block = z[z.group.eq(group)]
            old = cached[(cached.score_threshold.eq(threshold)) & cached.group.eq(group)]
            checks[f"{threshold}_{group}_count_and_KR_DE_match_saved_profiles"] = bool(
                len(block) == old.n_sequences.iloc[0] and all(np.isclose(block[f].mean(), old[old.feature.eq(f)]["mean"].iloc[0]) for f in ["KR_fraction", "DE_fraction"]))
        # For each anchor, reference probability of another sequence from the
        # same source/exact length sharing its group, excluding the anchor.
        total = z.groupby(["source", "length"]).group.transform("size")
        same = z.groupby(["source", "length", "group"]).group.transform("size")
        reference = (same - 1) / (total - 1).replace(0, np.nan)
        left = pairs.sequence_from
        right = pairs.sequence_to
        observed = z.group.reindex(left).to_numpy() == z.group.reindex(right).to_numpy()
        expected = (reference.reindex(left).to_numpy() + reference.reindex(right).to_numpy()) / 2
        assert np.isfinite(expected).all()
        neighbors.append({"threshold": threshold, "n_pairs": len(pairs), "same_group_n": int(observed.sum()),
                          "observed_fraction": observed.mean(), "source_length_reference_fraction": expected.mean(),
                          "same_vote_fraction": (votes.reindex(left).to_numpy() == votes.reindex(right).to_numpy()).mean(),
                          "identical_four_decisions_fraction": np.all(bits.reindex(left).to_numpy() == bits.reindex(right).to_numpy(), axis=1).mean()})
        # Explicit comparator lists check the probability formula independently.
        sample = np.linspace(0, len(pairs) - 1, 24, dtype=int)
        explicit = []
        for anchor in left.iloc[sample]:
            row = z.loc[anchor]
            candidates = z[(z.source == row.source) & (z.length == row.length)].drop(index=anchor)
            explicit.append(np.isclose(reference.loc[anchor], candidates.group.eq(row.group).mean()))
        checks[f"{threshold}_24_explicit_anchor_references_match"] = all(explicit)
        checks[f"{threshold}_votes_partition_full_cohort"] = sum(r["n"] for r in votes_rows if r["threshold"] == threshold) == len(z)
        checks[f"{threshold}_AMP_CLIP_decisions_unchanged"] = bool(bits[CLASS["AMP-CLIP"]].eq(scores[CLASS["AMP-CLIP"]].eq(1)).all())
        if previous_votes is not None:
            checks[f"{threshold}_higher_rule_cannot_increase_votes"] = bool(votes.le(previous_votes).all())
        previous_votes = votes
        for before in GROUPS:
            for after in GROUPS:
                transitions.append({"threshold": threshold, "baseline_group": before, "new_group": after,
                                    "n": int((frame.group.eq(before) & z.group.eq(after)).sum())})
    votes = pd.DataFrame(votes_rows)
    matched = pd.DataFrame(contrasts)
    neighbor = pd.DataFrame(neighbors)
    checks["baseline_neighbors_match_saved_result"] = int(neighbor[neighbor.threshold.eq(.5)].same_group_n.iloc[0]) == 32423
    checks["all_matched_length_differences_zero"] = bool(matched[matched.feature.eq("length")].right_minus_left.abs().lt(1e-12).all())
    assert all(checks.values()), checks
    for name, data in [("vote_counts", votes), ("matched_contrasts", matched), ("neighbor_similarity", neighbor), ("group_transitions", pd.DataFrame(transitions))]:
        data.to_csv(OUT / (name + ".csv"), index=False)
    pd.DataFrame([{"check": name, "passed": value} for name, value in checks.items()]).to_csv(OUT / "checks.csv", index=False)
    plot(votes, matched, neighbor)
    manifest = {"thresholds": RULES, "positive_rule": "score >= threshold; negative score < threshold",
                "n_complete": len(scores), "n_pairs": len(pairs), "AMP_CLIP_fixed_binary": True,
                "scope": "Sensitivity to two additional score rules; no calibrated probability, causal, or accuracy interpretation. Matched contrasts use equal sequence weights within source/exact-length cells; pairs remain dependent.",
                "checks": checks, "all_checks_passed": all(checks.values()),
                "input_hashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw_path, frame_path, pair_path]}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(neighbor.to_string(index=False))
    print(matched[matched.left.eq("All negative") & matched.right.eq("All positive") & matched.feature.isin(["KR_fraction", "DE_fraction"])][["threshold", "feature", "right_minus_left"]].to_string(index=False))
    print(f"All {len(checks)} checks passed.")


def plot(votes, matched, neighbor):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.7), layout="constrained")
    x = np.arange(3)
    bottom = np.zeros(3)
    for group, vote_set, color in [("All negative", [0], "#8a9da8"), ("Mixed", [1, 2, 3], "#bd6a34"), ("All positive", [4], "#287a65")]:
        values = votes[votes.positive_votes.isin(vote_set)].groupby("threshold").fraction_of_all.sum().reindex(RULES).to_numpy() * 100
        axes[0].bar(x, values, bottom=bottom, label=group, color=color, width=.6)
        for i, value in enumerate(values):
            axes[0].text(i, bottom[i] + value / 2, f"{value:.1f}%", ha="center", va="center", color="white", fontsize=9)
        bottom += values
    axes[0].set(ylabel="% of all 195,620 sequences", ylim=(0, 100), title="A. Decision groups")
    axes[0].legend(loc="upper center", bbox_to_anchor=(.5, -.17), fontsize=9, ncol=2)
    for feature, label, shift, color in [("KR_fraction", "K/R", -.17, "#225b78"), ("DE_fraction", "D/E", .17, "#bd6a34")]:
        rows = matched[matched.left.eq("All negative") & matched.right.eq("All positive") & matched.feature.eq(feature)].set_index("threshold")
        values = rows.reindex(RULES).right_minus_left.to_numpy() * 100
        axes[1].bar(x + shift, values, width=.30, label=label, color=color)
        for i, value in enumerate(values):
            axes[1].text(i + shift, value + (.5 if value > 0 else -.5), f"{value:.1f}", ha="center", va="bottom" if value > 0 else "top", fontsize=9)
    axes[1].axhline(0, color="#8a9da8", lw=.7)
    axes[1].set(ylabel="All positive − all negative (percentage points)", title="B. Within source and exact length", ylim=(-12, 20))
    axes[1].legend(loc="upper center", bbox_to_anchor=(.5, -.17), ncol=2)
    for feature, label, color in [("observed_fraction", "One-residue neighbors", "#225b78"), ("source_length_reference_fraction", "Source/length reference", "#bd6a34")]:
        values = neighbor.set_index("threshold").reindex(RULES)[feature].to_numpy() * 100
        axes[2].plot(x, values, "o-", color=color, label=label)
        for i, value in enumerate(values):
            axes[2].text(i, value + 1.5, f"{value:.1f}", ha="center", fontsize=9)
    axes[2].set(ylabel="% sharing the same decision group", ylim=(35, 95), title="C. Same 41,253 neighbor pairs")
    axes[2].legend(loc="upper center", bbox_to_anchor=(.5, -.17), fontsize=9)
    for ax in axes:
        ax.set(xticks=x, xticklabels=["0.3", "0.5", "0.7"], xlabel="Positive-score threshold")
    fig.suptitle("Do the findings depend on the positive/negative score rule?", fontsize=15, weight="bold")
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(OUT / ("classifier_rule_sensitivity." + ext), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
