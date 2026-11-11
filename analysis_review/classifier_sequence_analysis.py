"""Compare classifier consensus groups through composition and observed neighbors.

Uses saved predictions and the existing exact one-substitution pair catalogue.
No model fitting, new inference, downloads, or biological accuracy evaluation.
Run: python analysis_review/classifier_sequence_analysis.py
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/classifier_sequences"
CLASS = {"MACREL": "macrel_general_AMP", "AI4AMP": "ai4amp_general_AMP",
         "AMP-CLIP": "amp_clip_general_AMP", "AMPlify": "amplify_general_AMP"}
GROUPS = ["All negative", "Mixed", "All positive"]
AA = list("ACDEFGHIKLMNPQRSTVWY")
FEATURES = ["length", "KR_fraction", "DE_fraction", "charge_proxy", "hydrophobic_fraction",
            "contains_C", "entropy_bits"]


def group_labels(votes):
    return np.select([votes.eq(0), votes.eq(4)], [GROUPS[0], GROUPS[2]], default=GROUPS[1])


def save(data, name):
    pd.DataFrame(data).to_csv(OUT / (name + ".csv"), index=False)


def profile_rows(frame, by, population, weights=None):
    rows = []
    for label, block in frame.groupby(by, observed=True):
        w = np.ones(len(block)) if weights is None else weights.reindex(block.index).to_numpy()
        for feature in FEATURES:
            rows.append({"population": population, "grouping": by, "group": label, "feature": feature,
                         "n_sequences": len(block), "weight_sum": w.sum(),
                         "mean": np.average(block[feature], weights=w),
                         "median": block[feature].median(), "q25": block[feature].quantile(.25),
                         "q75": block[feature].quantile(.75),
                         "quantiles_weighting": "unweighted sequences"})
    return rows


def matched_contrasts(frame, weights, population):
    """Shared source/exact-length cells; min group mass weights each cell."""
    rows = []
    for left, right in combinations(GROUPS, 2):
        z = frame[frame.group.isin([left, right])].copy()
        z["mass"] = weights.reindex(z.index)
        keys = ["source", "length", "group"]
        mass = z.groupby(keys).mass.sum().unstack("group").reindex(columns=[left, right]).dropna()
        common = mass.index
        cell_weight = mass.min(axis=1)
        support = pd.MultiIndex.from_frame(z[["source", "length"]]).isin(common)
        for feature in FEATURES:
            weighted = z[feature] * z.mass
            means = weighted.groupby([z.source, z.length, z.group]).sum().unstack("group")
            means = means.reindex(index=common, columns=[left, right]) / mass
            delta = means[right] - means[left]
            rows.append({"population": population, "left": left, "right": right, "feature": feature,
                         "right_minus_left": np.average(delta, weights=cell_weight),
                         "unadjusted_right_minus_left": np.average(z.loc[z.group.eq(right), feature], weights=z.loc[z.group.eq(right), "mass"]) -
                             np.average(z.loc[z.group.eq(left), feature], weights=z.loc[z.group.eq(left), "mass"]),
                         "common_cells": len(common), "sum_minimum_cell_mass": cell_weight.sum(),
                         "left_n_in_common_cells": int((support & z.group.eq(left)).sum()),
                         "right_n_in_common_cells": int((support & z.group.eq(right)).sum()),
                         "left_total_n": int(z.group.eq(left).sum()), "right_total_n": int(z.group.eq(right).sum())})
    return rows


def neighbor_analysis(frame, scores, pair_path):
    pairs = pd.read_csv(pair_path, usecols=["sequence_from", "sequence_to", "source", "component"])
    pairs = pairs[pairs.sequence_from.isin(frame.index) & pairs.sequence_to.isin(frame.index)].copy()
    a = frame.loc[pairs.sequence_from].reset_index(drop=True)
    b = frame.loc[pairs.sequence_to].reset_index(drop=True)
    assert (a.length.to_numpy() == b.length.to_numpy()).all()
    assert (a.source.to_numpy() == b.source.to_numpy()).all()
    assert all(sum(x != y for x, y in zip(s, t)) == 1 for s, t in zip(pairs.sequence_from, pairs.sequence_to))
    pairs["group_from"] = a.group.to_numpy(); pairs["group_to"] = b.group.to_numpy()
    pairs["votes_from"] = a.votes.to_numpy(); pairs["votes_to"] = b.votes.to_numpy()
    bits_a = scores.loc[pairs.sequence_from].ge(.5).to_numpy()
    bits_b = scores.loc[pairs.sequence_to].ge(.5).to_numpy()
    pairs["classifier_decisions_changed"] = (bits_a != bits_b).sum(axis=1)
    same = a.group.to_numpy() == b.group.to_numpy()
    same_vote = a.votes.to_numpy() == b.votes.to_numpy()
    total = frame.groupby(["source", "length"]).size()
    counts = frame.groupby(["source", "length", "group"]).size()
    # Baseline excludes the anchor itself. A second reference excludes its entire
    # observed component; unavailable outside-component comparisons stay missing.
    mapped = frame[frame.component.notna()]
    ct = mapped.groupby(["source", "length", "component"]).size()
    cg = mapped.groupby(["source", "length", "component", "group"]).size()

    def reference(anchor, outside=False):
        key = pd.MultiIndex.from_frame(anchor[["source", "length"]])
        gkey = pd.MultiIndex.from_frame(anchor[["source", "length", "group"]])
        den = total.reindex(key).to_numpy(dtype=float)
        num = counts.reindex(gkey).to_numpy(dtype=float)
        if outside:
            ckey = pd.MultiIndex.from_frame(anchor[["source", "length", "component"]])
            cgkey = pd.MultiIndex.from_frame(anchor[["source", "length", "component", "group"]])
            den -= ct.reindex(ckey).to_numpy(dtype=float)
            num -= cg.reindex(cgkey).to_numpy(dtype=float)
        else:
            den -= 1; num -= 1
        return np.divide(num, den, out=np.full(len(anchor), np.nan), where=den > 0)

    ordinary = (reference(a) + reference(b)) / 2
    outside = (reference(a, True) + reference(b, True)) / 2
    pairs["same_group"] = same; pairs["same_vote_count"] = same_vote
    pairs["source_length_expected_same_group"] = ordinary
    pairs["outside_component_expected_same_group"] = outside
    pairs.to_csv(OUT / "one_substitution_pairs.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
    component_weight = 1 / pairs.component.map(pairs.component.value_counts()).to_numpy()
    hashed = pairs.apply(lambda row: hashlib.sha256((row.sequence_from + "|" + row.sequence_to).encode()).hexdigest(), axis=1)
    chosen = pairs.loc[hashed.sort_values().index].drop_duplicates("component").index
    populations = {"all_pairs": np.ones(len(pairs)), "equal_component_weight": component_weight,
                   "one_pair_per_component": pairs.index.isin(chosen).astype(float)}
    summary = []
    for pop, w in populations.items():
        for name, obs, ref in [("same_three_way_group", same, ordinary),
                               ("same_group_vs_outside_component", same, outside),
                               ("same_vote_count", same_vote, np.full(len(pairs), np.nan)),
                               ("identical_four_model_pattern", pairs.classifier_decisions_changed.eq(0).to_numpy(), np.full(len(pairs), np.nan))]:
            keep = w > 0
            if name == "same_group_vs_outside_component": keep &= np.isfinite(ref)
            summary.append({"population": pop, "metric": name, "n_pairs": int(keep.sum()),
                            "weight_sum": w[keep].sum(), "observed_fraction": np.average(obs[keep], weights=w[keep]),
                            "reference_fraction": np.average(ref[keep], weights=w[keep]) if np.isfinite(ref[keep]).all() else np.nan})
    save(summary, "neighbor_similarity")
    transitions = []
    anchors = pd.concat([pd.DataFrame({"anchor": a.group, "neighbor": b.group, "reference": reference(a)}),
                         pd.DataFrame({"anchor": b.group, "neighbor": a.group, "reference": reference(b)})], ignore_index=True)
    for label, block in anchors.groupby("anchor"):
        for neighbor in GROUPS:
            transitions.append({"anchor_group": label, "neighbor_group": neighbor, "directed_pair_n": len(block),
                                "count": int(block.neighbor.eq(neighbor).sum()), "fraction": block.neighbor.eq(neighbor).mean(),
                                "expected_same_group_source_length": block.reference.mean()})
    save(transitions, "neighbor_group_transitions")
    source_rows = []
    for source, ix in pairs.groupby("source").groups.items():
        source_rows.append({"source": source, "n_pairs": len(ix), "same_group_fraction": same[ix].mean(),
                            "source_length_reference_fraction": ordinary[ix].mean()})
    save(source_rows, "neighbor_by_source")
    save(pairs.groupby("classifier_decisions_changed").size().rename("n_pairs").reset_index(), "neighbor_decision_changes")
    return pairs


def main():
    start = time.monotonic(); OUT.mkdir(parents=True, exist_ok=True)
    raw_path = ROOT / "df_damp_195625x57.tsv"
    source_path = ROOT / "scratchpad/source_sequences.tsv"
    component_path = ROOT / "analysis_review/outputs/sanity_checks/broader_sequence_components.csv"
    pair_path = ROOT / "analysis_review/outputs/sequence_disagreement/single_substitution_pairs.csv.gz"
    scores = pd.read_csv(raw_path, sep="\t", usecols=["sequence", *CLASS.values()]).set_index("sequence", verify_integrity=True).dropna().sort_index()
    assert len(scores) == 195620 and scores.ge(0).all().all() and scores.le(1).all().all()
    source = pd.read_csv(source_path, sep="\t").set_index("sequence", verify_integrity=True).source
    components = pd.read_csv(component_path).set_index("sequence", verify_integrity=True).component
    seq = scores.index.to_series(); length = seq.str.len()
    residues = pd.DataFrame({aa: seq.str.count(aa) / length for aa in AA})
    frame = pd.DataFrame({"source": source.reindex(seq.index), "component": components.reindex(seq.index),
                          "length": length, "KR_fraction": residues.K + residues.R,
                          "DE_fraction": residues.D + residues.E,
                          "hydrophobic_fraction": residues[list("AILMFWVY")].sum(axis=1),
                          "contains_C": seq.str.contains("C").astype(float),
                          "entropy_bits": -(residues * np.log2(residues.where(residues.gt(0), 1))).sum(axis=1)})
    frame["charge_proxy"] = frame.KR_fraction - frame.DE_fraction
    frame["votes"] = scores.ge(.5).sum(axis=1); frame["group"] = group_labels(frame.votes)
    assert frame.source.notna().all() and np.allclose(residues.sum(axis=1), 1)
    frame.to_csv(OUT / "sequence_groups.csv.gz", index_label="sequence", compression={"method": "gzip", "mtime": 0})
    counts = frame.votes.value_counts().sort_index()
    mixed = int(counts.loc[1:3].sum())
    save([{"positive_votes": int(k), "n": int(n), "fraction_of_all": n / len(frame),
           "fraction_of_mixed": n / mixed if 0 < k < 4 else np.nan} for k, n in counts.items()], "vote_counts")
    profiles = profile_rows(frame, "group", "all_complete") + profile_rows(frame, "votes", "all_complete")
    mapped = frame[frame.component.notna()]
    weights = 1 / mapped.groupby(["source", "component"]).component.transform("size")
    profiles += profile_rows(mapped, "group", "mapped_equal_source_component", weights)
    save(profiles, "sequence_profiles")
    contrasts = matched_contrasts(frame, pd.Series(1., index=frame.index), "all_complete")
    contrasts += matched_contrasts(mapped, weights, "mapped_equal_source_component")
    save(contrasts, "source_length_matched_contrasts")
    per_source = []
    for label, block in frame.groupby("source"):
        per_source += [{"source": label, **r} for r in matched_contrasts(block, pd.Series(1., index=block.index), "within_source")]
    save(per_source, "within_source_length_contrasts")
    composition = []
    for label, ix in frame.groupby("group").groups.items():
        for aa in AA:
            composition.append({"group": label, "amino_acid": aa, "n_sequences": len(ix),
                                "mean_per_sequence_fraction": residues.loc[ix, aa].mean(),
                                "pooled_residue_fraction": (residues.loc[ix, aa] * length.loc[ix]).sum() / length.loc[ix].sum()})
    save(composition, "residue_composition")
    source_profiles = frame.groupby(["source", "group"])[FEATURES].mean().reset_index()
    source_profiles["n_sequences"] = frame.groupby(["source", "group"]).size().to_numpy()
    save(source_profiles, "source_profiles")
    sensitivity = []
    for threshold in [.3, .5, .7]:
        alt = frame.copy(); alt["group"] = group_labels(scores.ge(threshold).sum(axis=1))
        sensitivity += [{"score_threshold": threshold, **row} for row in profile_rows(alt, "group", "all_complete")]
    save(sensitivity, "score_rule_sensitivity")
    pairs = neighbor_analysis(frame, scores, pair_path)
    expected = pd.read_csv(ROOT / "analysis_review/outputs/screening_decision_trends/classifier_vote_patterns.csv", dtype={"pattern": str})
    expected = expected[expected.population.eq("all_classifier_complete")]
    check_counts = expected.groupby(expected.pattern.str.zfill(4).map(lambda x: x.count("1")))["count"].sum()
    checks = {"195620_classifier_complete_sequences": len(frame) == 195620,
              "vote_counts_match_independently_saved_16_patterns": counts.equals(check_counts.rename_axis(counts.index.name)),
              "three_way_partition_matches_existing_counts": frame.group.value_counts().to_dict() == {"Mixed": 128097, "All positive": 42645, "All negative": 24878},
              "residue_fractions_sum_to_one": bool(np.allclose(residues.sum(axis=1), 1)),
              "exact_length_matched_contrasts_are_zero": all(abs(r["right_minus_left"]) < 1e-10 for r in contrasts if r["feature"] == "length"),
              "one_substitution_catalogue_retained": len(pairs) == 41253,
              "pairs_same_source_and_single_substitution": True,
              "no_duplicate_unordered_pairs": not pairs[["sequence_from", "sequence_to"]].duplicated().any()}
    save([{"check": k, "passed": bool(v)} for k, v in checks.items()], "checks")
    assert all(checks.values()), checks
    manifest = {"scope": "Exploratory classifier-group composition and sequence-neighbor analysis; no biological accuracy claim",
                "date": "2026-11-10", "n_complete": len(frame), "n_component_mapped": len(mapped),
                "software_versions": {"numpy": np.__version__, "pandas": pd.__version__},
                "n_component_unmapped": len(frame) - len(mapped), "models": CLASS, "positive_rule": "score >= 0.5; AMP-CLIP is already binary",
                "neighbor_scope": "Previously enumerated same-source one-substitution pairs in the MIC-complete cohort, not exhaustive classifier-wide neighbors",
                "matched_estimate": "Within source/exact-length mean differences, weighted by smaller group mass; no matching on model outputs",
                "inputs_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw_path, source_path, component_path, pair_path]},
                "checks": checks, "elapsed_seconds": time.monotonic() - start}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Classifier analysis: {len(frame):,} sequences, {len(pairs):,} local pairs; all {len(checks)} checks passed.")


if __name__ == "__main__":
    main()
