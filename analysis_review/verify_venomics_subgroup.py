"""Separate raw-data arithmetic checks for the Venomics report.

Run after venomics_subgroup.py. Uses NumPy selection and explicit stratum loops,
instead of the analysis script's pandas selection and grouped transforms.
"""
from pathlib import Path
from collections import Counter
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/venomics"
MIC = ["apex_1_1_E._coli_ATCC_11775", "mbc_attention_E._coli_unspecified_strain", "ampredictor_E._coli_unspecified_strain"]
CLASS = ["macrel_general_AMP", "ai4amp_general_AMP", "amp_clip_general_AMP", "amplify_general_AMP"]


def main():
    raw_path = ROOT / "df_damp_195625x57.tsv"
    raw = pd.read_csv(raw_path, sep="\t").set_index("sequence", verify_integrity=True)
    raw = raw.dropna(subset=MIC).sort_index()
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence").source.reindex(raw.index)
    component = pd.read_csv(ROOT / "analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence").component.reindex(raw.index)
    context = pd.read_csv(OUT / "source_context.csv").set_index("source")
    annotations = pd.read_csv(OUT / "venomics_sequence_annotations.csv").set_index("sequence")
    v = source.eq("Venomics")
    votes = np.sum(raw[CLASS].to_numpy() >= .5, axis=1)
    z = pd.DataFrame({"source": source, "component": component, "length": raw.index.str.len(),
                      "contains_C": raw.index.str.contains("C").astype(float),
                      "KR_fraction": [(s.count("K") + s.count("R")) / len(s) for s in raw.index],
                      "all_positive": (votes == 4).astype(float),
                      "APEX_pass16": raw[MIC[0]].le(16).astype(float),
                      "MBC_pass16": raw[MIC[1]].le(16).astype(float),
                      "AMP_pass16": raw[MIC[2]].le(16).astype(float)})
    checks = {"primary_cohort_and_Venomics_count": len(raw) == 194947 and int(v.sum()) == 4618,
              "source_counts_match_context": context.n.to_dict() == source.value_counts().to_dict(),
              "raw_votes_2536_positive_2043_mixed_39_negative": [int(np.sum(votes[v] == 4)), int(np.sum((votes[v] > 0) & (votes[v] < 4))), int(np.sum(votes[v] == 0))] == [2536, 2043, 39],
              "sequence_features_and_components_match_annotations": bool(np.allclose(z.loc[v, ["KR_fraction", "contains_C"]], annotations.loc[z.index[v], ["KR_fraction", "contains_C"]]) and np.array_equal(z.loc[v, "component"].to_numpy(), annotations.loc[z.index[v], "component"].to_numpy()))}
    group_sizes = Counter(zip(source, component))
    z["weight"] = [1 / group_sizes[(s, c)] for s, c in zip(source, component)]
    for s in context.index:
        mask = source.eq(s)
        c = Counter(component[mask])
        checks[f"{s}_groups_and_top10_recomputed"] = len(c) == context.loc[s, "n_groups"] and np.isclose(sum(sorted(c.values(), reverse=True)[:10]) / mask.sum(), context.loc[s, "top10_groups_fraction"])
        for metric in ["all_positive", "APEX_pass16", "MBC_pass16", "AMP_pass16"]:
            checks[f"{s}_{metric}_raw_and_group_weighted"] = bool(np.isclose(z.loc[mask, metric].mean(), context.loc[s, metric]) and np.isclose(np.average(z.loc[mask, metric], weights=z.loc[mask, "weight"]), context.loc[s, metric + "_equal_group"]))
    # Exact set membership from independently sorted raw MIC arrays.
    ven = raw.loc[v, MIC]
    chosen = np.zeros((len(ven), 3), dtype=bool)
    for j in range(3): chosen[np.argsort(ven.iloc[:, j].to_numpy(), kind="stable")[:231], j] = True
    patterns = np.bincount(chosen @ np.array([4, 2, 1]), minlength=8)
    expected = pd.read_csv(OUT / "venomics_exact_overlaps.csv", dtype={"pattern": str})
    cells = expected[expected.rule.eq("top_5%")]
    checks["all_top5_exclusive_cells_match_NumPy_sort"] = all(patterns[int(r.pattern, 2)] == r.n for r in cells.itertuples())
    cutoff = ven.to_numpy() <= 16
    patterns = np.bincount(cutoff @ np.array([4, 2, 1]), minlength=8)
    checks["all_cutoff_cells_match_raw_bit_counts"] = all(patterns[int(r.pattern, 2)] == r.n for r in expected[expected.rule.eq("MIC_le_16")].itertuples())
    # Every matched cell is explicitly visited, using raw sequence features.
    saved = pd.read_csv(OUT / "controlled_comparisons.csv")
    for name, add_kr, weighted in [("length_C", False, False), ("length_C_KR", True, False), ("equal_group_length_C_KR", True, True)]:
        keys = ["length", "contains_C"]
        z["KR_bin"] = np.floor(z.KR_fraction * 10 + 1e-10).astype(int)
        if add_kr: keys += ["KR_bin"]
        a, b, mass = [], [], []
        for _, block in z.groupby(keys):
            left = block[block.source.eq("Venomics")]
            right = block[block.source.ne("Venomics")]
            if left.empty or right.empty: continue
            metrics = ["all_positive", "APEX_pass16", "MBC_pass16", "AMP_pass16"]
            wl = left.weight if weighted else np.ones(len(left))
            wr = right.weight if weighted else np.ones(len(right))
            a.append(np.average(left[metrics], weights=wl, axis=0)); b.append(np.average(right[metrics], weights=wr, axis=0)); mass.append(np.sum(wl))
        delta = 100 * (np.average(a, weights=mass, axis=0) - np.average(b, weights=mass, axis=0))
        target = saved[saved.comparison.eq(name)].set_index("metric").loc[metrics]
        checks[name + "_all_four_outcomes_verified_by_explicit_cells"] = bool(np.allclose(delta, target.difference_percentage_points, atol=1e-10, rtol=0) and np.allclose(sum(mass), target.Venomics_mass_in_common_cells))
    old = pd.read_csv(ROOT / "analysis_review/outputs/shortlist_overlap/spearman_correlations.csv")
    new = pd.read_csv(OUT / "source_correlations.csv")
    joined = new.merge(old, left_on=["source", "left", "right"], right_on=["population", "left", "right"], suffixes=("_new", "_old"))
    checks["all_24_source_correlations_match_prior_review"] = len(joined) == 24 and bool(np.allclose(joined.spearman_new, joined.spearman_old, atol=1e-12))
    assert all(checks.values()), {k: bool(v) for k, v in checks.items() if not v}
    result = {"all_checks_passed": True, "n_checks": len(checks), "checks": {k: bool(v) for k, v in checks.items()},
              "raw_input_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
              "scope": "Separate arithmetic verification using the same saved predictions; no independent biological evidence."}
    (OUT / "independent_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Venomics verification: all {len(checks)} checks passed.")


if __name__ == "__main__":
    main()
