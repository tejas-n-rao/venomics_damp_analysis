"""Verify the source analysis through separate raw-array and explicit-cell checks."""
from pathlib import Path
from collections import Counter
import hashlib
import json

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/source_comparisons"
MIC = ["apex_1_1_E._coli_ATCC_11775", "mbc_attention_E._coli_unspecified_strain", "ampredictor_E._coli_unspecified_strain"]
NAMES = ["APEX 1.1", "MBC-Attention", "AMPredictor"]
CLASS = ["macrel_general_AMP", "ai4amp_general_AMP", "amp_clip_general_AMP", "amplify_general_AMP"]


def main():
    path = ROOT / "df_damp_195625x57.tsv"
    raw = pd.read_csv(path, sep="\t").set_index("sequence", verify_integrity=True).sort_index()
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence").source
    group = pd.read_csv(ROOT / "analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence").component
    primary = raw.dropna(subset=MIC)
    short = pd.read_csv(OUT / "shortlist_comparisons.csv")
    checks = {}
    for label in sorted(source.unique()):
        index = primary.index[source.reindex(primary.index).eq(label)]
        data = primary.loc[index, MIC].to_numpy()
        for budget in [.01, .02, .05, .1]:
            k = int(np.ceil(budget * len(data)))
            chosen = [set(np.argsort(data[:, j], kind="stable")[:k]) for j in range(3)]
            expected = short[short.source.eq(label) & short.population.eq("every_sequence") & short.fraction.eq(budget)]
            checks[f"{label}_budget_{budget}_NumPy_intersections"] = all(len(chosen[NAMES.index(r.left)] & chosen[NAMES.index(r.right)]) == r.shared_n for r in expected.itertuples())
        # Independently pick the minimum hash in each group, without sorting a frame.
        choices = {}
        for sequence, component in group.reindex(index).items():
            value = hashlib.sha256(f"0|{sequence}".encode()).hexdigest()
            if component not in choices or value < choices[component][0]: choices[component] = (value, sequence)
        selected = primary.loc[[item[1] for item in choices.values()], MIC]
        k = int(np.ceil(len(selected) * .05))
        model_sets = [set(selected.index[np.argsort(selected.iloc[:, j].to_numpy(), kind="stable")[:k]]) for j in range(3)]
        expected = short[short.source.eq(label) & short.population.eq("representative_00")]
        checks[f"{label}_minimum_hash_representatives"] = all(len(model_sets[NAMES.index(r.left)] & model_sets[NAMES.index(r.right)]) == r.shared_n for r in expected.itertuples())
    decisions = pd.read_csv(OUT / "matched_target_decisions.csv")
    specs = {"E. coli": ["deep_amp_cnn_neg_Gram-negative_unspecified_species_and_strain", "deep_amp_lstm_neg_Gram-negative_unspecified_species_and_strain"],
             "B. subtilis": ["deep_amp_cnn_pos_Gram-positive_unspecified_species_and_strain", "deep_amp_lstm_pos_Gram-positive_unspecified_species_and_strain"],
             "APEX E. coli versions": ["apex_1_0_E._coli_ATCC11775", MIC[0]]}
    for row in decisions.itertuples():
        p = raw.loc[source.index[source.eq(row.source)], specs[row.target]].dropna()
        array = p.to_numpy() <= row.cutoff_uM
        cells = np.bincount(array @ np.array([2, 1]), minlength=4)
        checks[f"{row.source}_{row.target}_{row.cutoff_uM}_four_cells"] = (cells[0] == row.both_fail_n and cells[3] == row.both_pass_n and cells[1] + cells[2] == row.disagree_n and len(p) == row.n_paired)
    # Explicit length cells check both raw and group-weighted KR/DE contrasts.
    z = pd.DataFrame({"source": source.reindex(primary.index), "length": primary.index.str.len(), "component": group.reindex(primary.index),
                      "KR_fraction": [(s.count("K") + s.count("R")) / len(s) for s in primary.index],
                      "DE_fraction": [(s.count("D") + s.count("E")) / len(s) for s in primary.index]})
    votes = (primary[CLASS].to_numpy() >= .5).sum(axis=1)
    z["label"] = np.where(votes == 4, "positive", np.where(votes == 0, "negative", "mixed"))
    profiles = pd.read_csv(OUT / "classifier_length_matched_profiles.csv")
    for label, f in z.groupby("source"):
        sizes = Counter(f.component)
        for population in ["every_sequence", "equal_group"]:
            weight = pd.Series(1., index=f.index) if population == "every_sequence" else f.component.map(lambda g: 1 / sizes[g])
            differences, mass = [], []
            for _, block in f.groupby("length"):
                a = block[block.label.eq("negative")]; b = block[block.label.eq("positive")]
                if a.empty or b.empty: continue
                wa = weight.loc[a.index]; wb = weight.loc[b.index]
                differences.append(np.average(b[["KR_fraction", "DE_fraction"]], weights=wb, axis=0) - np.average(a[["KR_fraction", "DE_fraction"]], weights=wa, axis=0))
                mass.append(min(wa.sum(), wb.sum()))
            actual = np.average(differences, weights=mass, axis=0)
            target = profiles[profiles.source.eq(label) & profiles.population.eq(population)].set_index("feature").loc[["KR_fraction", "DE_fraction"]]
            checks[f"{label}_{population}_explicit_KR_DE_contrast"] = bool(np.allclose(actual, target.right_minus_left, atol=1e-12, rtol=0) and np.allclose(sum(mass), target.sum_minimum_cell_mass))
    # Ratio >=2 is checked separately from the analysis script's log2-delta >=1.
    pairs = pd.read_csv(ROOT / "analysis_review/outputs/screening_decision_trends/DeepAMP_composition_controlled_pairs.csv")
    order = pd.read_csv(OUT / "sequence_order_responses.csv")
    for row in order.itertuples():
        p = pairs[pairs.source.eq(row.source)]
        va = raw.loc[p.sequence_left, specs[row.target]].to_numpy(); vb = raw.loc[p.sequence_right, specs[row.target]].to_numpy()
        if row.population == "no_group_reuse": keep = p.greedy_component_representative.to_numpy(dtype=bool)
        elif row.population == "all_predictions_le128": keep = np.all((va <= 128) & (vb <= 128), axis=1)
        else: keep = np.ones(len(p), bool)
        j = 0 if row.model == "CNN" else 1
        ratio = np.maximum(va[keep, j], vb[keep, j]) / np.minimum(va[keep, j], vb[keep, j])
        checks[f"{row.source}_{row.target}_{row.population}_{row.model}_ratio_check"] = int(keep.sum()) == row.n_pairs and int(np.sum(ratio >= 2)) == row.twofold_n
    chosen = pairs[pairs.greedy_component_representative]
    used = set(); no_reuse = True
    for a, b in zip(chosen.sequence_left, chosen.sequence_right):
        ids = {group[a], group[b]}
        no_reuse &= not bool(ids & used); used.update(ids)
    checks["2197_retained_order_pairs_have_no_group_reuse"] = no_reuse and len(chosen) == 2197
    previous = pd.read_csv(ROOT / "analysis_review/outputs/classifier_sequences/neighbor_by_source.csv").set_index("source")
    current = pd.read_csv(OUT / "neighbor_agreement.csv").query('population == "every_pair"').set_index("source")
    checks["all_eight_observed_neighbor_rates_reproduce"] = bool(np.allclose(previous.same_group_fraction, current.same_group_fraction.reindex(previous.index), atol=1e-12))
    assert all(checks.values()), {k: bool(v) for k, v in checks.items() if not v}
    result = {"all_checks_passed": True, "n_checks": len(checks), "checks": {k: bool(v) for k, v in checks.items()},
              "raw_input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "scope": "Separate arithmetic implementations on the same predictions, not independent biological validation."}
    (OUT / "independent_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"Source-stratified verification: all {len(checks)} checks passed.")


if __name__ == "__main__":
    main()
