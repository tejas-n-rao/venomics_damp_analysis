"""DAMP screening decisions with concentration cutoffs and matched output targets.

Run with the existing environment: python analysis_review/screening_decision_trends.py
No model inference, input edits, calibration fit, or experimental accuracy claim.
Cutoffs and margins are exploratory operational rules, not clinical breakpoints.
Native Gram labels are retained: the paper reports E. coli and B. subtilis
training for the negative/positive heads, with synthetic log-MIC labels for
non-AMP examples. These labels do not establish whole-Gram-class activity.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import os
import tempfile
import time

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/screening_decision_trends"
OUT.mkdir(parents=True, exist_ok=True)
MIC = {"APEX 1.1": "apex_1_1_E._coli_ATCC_11775",
       "MBC-Attention": "mbc_attention_E._coli_unspecified_strain",
       "AMPredictor": "ampredictor_E._coli_unspecified_strain"}
CLASS = {"MACREL": "macrel_general_AMP", "AI4AMP": "ai4amp_general_AMP",
         "AMP-CLIP": "amp_clip_general_AMP", "AMPlify": "amplify_general_AMP"}
CUTOFFS = [4, 8, 16, 32, 64, 128]
CHECKS = []


def check(name, condition):
    CHECKS.append({"check": name, "passed": bool(condition)})


def save(data, name):
    pd.DataFrame(data).to_csv(OUT / (name + ".csv"), index=False)


def fraction(a, b):
    return a / b if b else np.nan


def decision_pair(a, b, cutoff):
    aa = a.le(cutoff); bb = b.le(cutoff)
    both = int((aa & bb).sum()); only_a = int((aa & ~bb).sum()); only_b = int((~aa & bb).sum())
    neither = int((~aa & ~bb).sum()); n = len(a)
    return {"n": n, "left_pass_n": int(aa.sum()), "right_pass_n": int(bb.sum()),
            "both_pass_n": both, "left_only_n": only_a, "right_only_n": only_b, "neither_pass_n": neither,
            "decision_disagreement_n": only_a + only_b,
            "decision_disagreement_fraction": (only_a + only_b) / n,
            "both_fraction_of_either_passing": fraction(both, both + only_a + only_b),
            "right_pass_fraction_of_left_pass": fraction(both, both + only_a),
            "left_pass_fraction_of_right_pass": fraction(both, both + only_b)}


def representative(index, components):
    eligible = index[components.reindex(index).notna()]
    order = sorted(eligible, key=lambda s: hashlib.sha256(s.encode()).hexdigest())
    ids = components.reindex(order)
    return pd.Index(order)[~ids.duplicated()]


def main():
    started = time.monotonic(); path = ROOT / "df_damp_195625x57.tsv"
    raw = pd.read_csv(path, sep="\t").set_index("sequence", verify_integrity=True)
    metadata = pd.DataFrame(json.loads((ROOT / "scratchpad/damp_metadata.json").read_text())).set_index("column")
    source = pd.read_csv(ROOT / "scratchpad/source_sequences.tsv", sep="\t").set_index("sequence", verify_integrity=True).source
    components = pd.read_csv(ROOT / "analysis_review/outputs/sanity_checks/broader_sequence_components.csv").set_index("sequence", verify_integrity=True).component
    length = raw.index.to_series().str.len(); cysteine = raw.index.to_series().str.contains("C")
    domain = length.between(8, 39) & ~cysteine
    mic = raw[list(MIC.values())].rename(columns={v: k for k, v in MIC.items()}).dropna().sort_index()
    populations = {"all_common": mic.index, "cysteine_free_length8_to39": mic.index[domain.reindex(mic.index)],
                   "one_component_representative": representative(mic.index, components)}
    check("three MIC outputs have common expected denominator", len(mic) == 194947)
    check("all common MICs finite and positive", np.isfinite(mic).all().all() and mic.gt(0).all().all())
    check("source and component map cover MIC cohort", source.reindex(mic.index).notna().all() and components.reindex(mic.index).notna().all())
    check("all MIC columns declared uM", metadata.loc[list(MIC.values()), "output_kind"].eq("MIC_uM").all())
    print("Analyzing fixed-concentration screening decisions", flush=True)
    counts = []; patterns = []; pairs = []; margins = []; by_source = []
    for pop, index in populations.items():
        x = mic.loc[index]
        previous = np.zeros(3, dtype=int)
        for cutoff in CUTOFFS:
            passed = x.le(cutoff); vote = passed.sum(axis=1)
            for model in x:
                counts.append({"population": pop, "cutoff_uM": cutoff, "model": model, "n": len(x),
                               "pass_n": int(passed[model].sum()), "pass_fraction": passed[model].mean(),
                               "sole_passing_model_n": int((passed[model] & vote.eq(1)).sum()),
                               "all_three_pass_n": int(vote.eq(3).sum()), "any_pass_n": int(vote.gt(0).sum())})
            codes = passed.astype(int).astype(str).agg("".join, axis=1).value_counts()
            for code in [f"{i:03b}" for i in range(8)]:
                patterns.append({"population": pop, "cutoff_uM": cutoff, "model_order": ";".join(x),
                                 "pattern": code, "n": len(x), "count": int(codes.get(code, 0)),
                                 "fraction": codes.get(code, 0) / len(x)})
            check(f"all threshold patterns reconcile {pop}/{cutoff}", codes.sum() == len(x))
            check(f"passing counts increase with cutoff {pop}/{cutoff}", (passed.sum().to_numpy() >= previous).all())
            previous = passed.sum().to_numpy()
            for left, right in combinations(x, 2):
                result = decision_pair(x[left], x[right], cutoff)
                pairs.append({"population": pop, "cutoff_uM": cutoff, "left": left, "right": right, **result})
                for margin in [1., np.sqrt(2), 2.]:
                    far = ((x[left].le(cutoff / margin) & x[right].gt(cutoff * margin)) |
                           (x[right].le(cutoff / margin) & x[left].gt(cutoff * margin)))
                    margins.append({"comparison_type": "approximate_species", "population": pop, "target": "E. coli",
                                    "left": left, "right": right, "cutoff_uM": cutoff, "margin_factor": margin,
                                    "n": len(x), "far_disagreement_n": int(far.sum()),
                                    "far_disagreement_fraction": far.mean(),
                                    "fraction_of_all_cutoff_disagreements": fraction(far.sum(), result["decision_disagreement_n"])})
                    if margin == 1: check(f"margin=1 reproduces disagreement {pop}/{cutoff}/{left}/{right}", far.sum() == result["decision_disagreement_n"])
            if pop == "all_common":
                for src, ix in source.reindex(index).groupby(source.reindex(index)).groups.items():
                    for model in x:
                        by_source.append({"source": src, "cutoff_uM": cutoff, "model": model, "n": len(ix),
                                          "pass_n": int(passed.loc[ix, model].sum()), "pass_fraction": passed.loc[ix, model].mean()})
    save(counts, "MIC_cutoff_counts"); save(patterns, "MIC_cutoff_vote_patterns")
    save(pairs, "MIC_cutoff_pairwise"); save(by_source, "MIC_cutoff_by_source")

    print("Analyzing target-matched comparisons", flush=True)
    matched_specs = []
    old = metadata[metadata.model.eq("apex_1_0") & ~metadata.alias_unverified]
    new = metadata[metadata.model.eq("apex_1_1") & ~metadata.alias_unverified]
    for col, row in new.iterrows():
        matches = old.index[old.strain_index.eq(row.strain_index)]
        if len(matches) == 1:
            matched_specs.append(("APEX_version_same_strain", row.target, matches[0], col))
    check("eight unambiguous same-strain APEX comparisons", len(matched_specs) == 8)
    for gram, target in [("neg", "Gram-negative"), ("pos", "Gram-positive")]:
        suffix = gram + "_" + target + "_unspecified_species_and_strain"
        matched_specs.append(("DeepAMP_same_Gram_group", target, "deep_amp_cnn_" + suffix, "deep_amp_lstm_" + suffix))
    matched = []; matched_sources = []; coverage = []; shifts = []; deep_features = []
    for kind, target, left, right in matched_specs:
        x = raw[[left, right]].dropna().sort_index()
        check("finite positive matched outputs " + target, np.isfinite(x).all().all() and x.gt(0).all().all())
        coverage.append({"comparison": kind, "target": target, "paired_n": len(x),
                         "paired_with_component_map_n": int(components.reindex(x.index).notna().sum())})
        paired_pops = {"all_paired": x.index, "cysteine_free_length8_to39": x.index[domain.reindex(x.index)],
                       "one_component_representative": representative(x.index, components)}
        for pop, ix in paired_pops.items():
            xx = x.loc[ix]; ratio = xx[right] / xx[left]
            shifts.append({"comparison": kind, "target": target, "population": pop, "left": left, "right": right,
                           "n": len(xx), "median_right_over_left_MIC": ratio.median(),
                           "fraction_right_MIC_lower": ratio.lt(1).mean(), "fraction_within_twofold": ratio.between(.5, 2).mean(),
                           "spearman": xx[left].corr(xx[right], method="spearman")})
            for cutoff in CUTOFFS:
                result = decision_pair(xx[left], xx[right], cutoff)
                matched.append({"comparison": kind, "target": target, "population": pop, "cutoff_uM": cutoff,
                                "left": left, "right": right, **result})
                for margin in [1., np.sqrt(2), 2.]:
                    far = ((xx[left].le(cutoff / margin) & xx[right].gt(cutoff * margin)) |
                           (xx[right].le(cutoff / margin) & xx[left].gt(cutoff * margin)))
                    margins.append({"comparison_type": kind, "population": pop, "target": target,
                                    "left": left, "right": right, "cutoff_uM": cutoff, "margin_factor": margin,
                                    "n": len(xx), "far_disagreement_n": int(far.sum()), "far_disagreement_fraction": far.mean(),
                                    "fraction_of_all_cutoff_disagreements": fraction(far.sum(), result["decision_disagreement_n"])})
            if pop == "all_paired":
                for src, idx in source.reindex(ix).groupby(source.reindex(ix)).groups.items():
                    matched_sources.append({"comparison": kind, "target": target, "source": src, "cutoff_uM": 16,
                                            **decision_pair(xx.loc[idx, left], xx.loc[idx, right], 16)})
        if kind == "DeepAMP_same_Gram_group":
            # Length/source controlled, descriptive sequence profile of direction
            # of disagreement at 16 uM. No comparison to strain-specific MICs.
            seq = x.index.to_series(); ll = seq.str.len()
            ff = pd.DataFrame({"source": source.reindex(x.index), "length": ll,
                               "KR_fraction": seq.str.count("[KR]") / ll,
                               "DE_fraction": seq.str.count("[DE]") / ll,
                               "hydrophobic_fraction": seq.str.count("[AILMFWVY]") / ll,
                               "contains_C": seq.str.contains("C").astype(float)})
            cnn_only = x[left].le(16) & x[right].gt(16)
            lstm_only = x[right].le(16) & x[left].gt(16)
            groups = {"CNN_only_passes": cnn_only, "LSTM_only_passes": lstm_only}
            for group, mask in groups.items():
                for feature in ["length", "KR_fraction", "DE_fraction", "hydrophobic_fraction", "contains_C"]:
                    deep_features.append({"target": target, "comparison": "unadjusted_profile", "group": group,
                                          "feature": feature, "n": int(mask.sum()), "mean": ff.loc[mask, feature].mean(),
                                          "median": ff.loc[mask, feature].median()})
            bb = ff.loc[cnn_only | lstm_only].copy()
            bb["group"] = np.where(cnn_only.loc[bb.index], "CNN_only", "LSTM_only")
            nn = bb.groupby(["source", "length", "group"]).size().unstack("group").dropna()
            mass = nn.min(axis=1)
            for feature in ["length", "KR_fraction", "DE_fraction", "hydrophobic_fraction", "contains_C"]:
                mm = bb.groupby(["source", "length", "group"])[feature].mean().unstack("group").loc[mass.index]
                delta = np.average(mm.CNN_only - mm.LSTM_only, weights=mass)
                deep_features.append({"target": target, "comparison": "same_source_exact_length", "group": "CNN_only_minus_LSTM_only",
                                      "feature": feature, "n": int(mass.sum()), "mean": delta, "median": np.nan,
                                      "common_strata": len(mass), "CNN_only_in_common_strata_n": int(nn.CNN_only.sum()),
                                      "LSTM_only_in_common_strata_n": int(nn.LSTM_only.sum())})
                if feature == "length": check("DeepAMP exact length matching " + target, np.isclose(delta, 0))
    save(matched, "matched_target_cutoff_decisions"); save(shifts, "matched_target_output_shifts")
    save(matched_sources, "matched_target_by_source"); save(coverage, "matched_target_coverage")
    save(deep_features, "DeepAMP_directional_sequence_profiles"); save(margins, "cutoff_margin_sensitivity")

    print("Analyzing classifier decisions and concentration cutoffs", flush=True)
    cp = raw[list(CLASS.values())].rename(columns={v: k for k, v in CLASS.items()}).dropna().sort_index()
    check("classifier exports in zero-one range", cp.ge(0).all().all() and cp.le(1).all().all())
    check("AMP-CLIP is exported as binary", set(cp["AMP-CLIP"].unique()) == {0., 1.})
    class_counts = []; class_patterns = []; class_pairs = []; dissent = []; bridge = []; class_source = []; class_references = []
    class_pops = {"all_classifier_complete": cp.index, "MIC_common": mic.index,
                  "cysteine_free_length8_to39": cp.index[domain.reindex(cp.index)],
                  "one_component_representative": representative(cp.index, components)}
    for pop, ix in class_pops.items():
        for threshold in [.3, .5, .7]:
            b = cp.loc[ix].ge(threshold); votes = b.sum(axis=1)
            # Expected unanimity under independent labels with the same model
            # positive rates. The stratified reference preserves source/length
            # rates; neither reference supplies ground truth or significance.
            pp = b.mean()
            pp_stratum = b.groupby([source.reindex(ix), b.index.str.len()]).transform("mean")
            class_references.append({"population": pop, "score_threshold": threshold, "n": len(b),
                                     "observed_unanimity_fraction": (votes.eq(0) | votes.eq(4)).mean(),
                                     "independent_marginal_expected_unanimity": pp.prod() + (1 - pp).prod(),
                                     "independent_within_source_length_expected_unanimity": (pp_stratum.prod(axis=1) + (1 - pp_stratum).prod(axis=1)).mean()})
            for model in b:
                class_counts.append({"population": pop, "score_threshold": threshold, "model": model, "n": len(b),
                                     "positive_n": int(b[model].sum()), "positive_fraction": b[model].mean(),
                                     "all_positive_n": int(votes.eq(4).sum()), "all_negative_n": int(votes.eq(0).sum()),
                                     "mixed_n": int(votes.between(1, 3).sum())})
                others = b.drop(columns=model).sum(axis=1)
                dissent.append({"population": pop, "score_threshold": threshold, "model": model,
                                "only_positive_n": int((b[model] & others.eq(0)).sum()),
                                "only_negative_n": int((~b[model] & others.eq(3)).sum())})
            codes = b.astype(int).astype(str).agg("".join, axis=1).value_counts()
            check(f"classifier patterns reconcile {pop}/{threshold}", int(codes.sum()) == len(b))
            if threshold == .5:
                for code in [f"{i:04b}" for i in range(16)]:
                    class_patterns.append({"population": pop, "model_order": ";".join(b), "pattern": code,
                                           "n": len(b), "count": int(codes.get(code, 0)), "fraction": codes.get(code, 0) / len(b)})
                for left, right in combinations(b, 2):
                    a = b[left]; c = b[right]
                    both = int((a & c).sum()); neither = int((~a & ~c).sum())
                    observed = a.eq(c).mean(); expected = a.mean() * c.mean() + (1 - a.mean()) * (1 - c.mean())
                    class_pairs.append({"population": pop, "left": left, "right": right, "n": len(b),
                                        "agreement_fraction": observed, "cohen_kappa": (observed - expected) / (1 - expected),
                                        "positive_agreement": fraction(2 * both, int(a.sum() + c.sum())),
                                        "negative_agreement": fraction(2 * neither, int((~a).sum() + (~c).sum()))})
                if pop == "all_classifier_complete":
                    for src, idx in source.reindex(ix).groupby(source.reindex(ix)).groups.items():
                        vv = votes.loc[idx]
                        class_source.append({"source": src, "n": len(idx), "all_positive_fraction": vv.eq(4).mean(),
                                             "all_negative_fraction": vv.eq(0).mean(), "mixed_fraction": vv.between(1, 3).mean()})
    votes = cp.loc[mic.index].ge(.5).sum(axis=1)
    for cutoff in CUTOFFS:
        for model in mic:
            selected = mic[model].le(cutoff)
            bridge.append({"cutoff_uM": cutoff, "MIC_model": model, "n": len(mic), "MIC_pass_n": int(selected.sum()),
                           "MIC_pass_all_classifiers_positive_n": int((selected & votes.eq(4)).sum()),
                           "MIC_pass_all_classifiers_negative_n": int((selected & votes.eq(0)).sum()),
                           "all_positive_fraction_among_MIC_pass": votes.loc[selected].eq(4).mean(),
                           "all_negative_fraction_among_MIC_pass": votes.loc[selected].eq(0).mean()})
    save(class_counts, "classifier_decision_counts"); save(class_patterns, "classifier_vote_patterns")
    save(class_pairs, "classifier_pairwise_agreement"); save(dissent, "classifier_single_dissent")
    save(class_source, "classifier_consensus_by_source"); save(bridge, "MIC_cutoff_classifier_bridge")
    save(class_references, "classifier_unanimity_independence_references")

    print("Checking composition-controlled DeepAMP order responses", flush=True)
    observed_pairs = pd.read_csv(ROOT / "analysis_review/outputs/additional_trends/exact_composition_two_member_pairs.csv")
    deep_columns = [c for c in raw if c.startswith("deep_amp_")]
    valid = (raw[deep_columns].reindex(observed_pairs.sequence_left).notna().all(axis=1).to_numpy() &
             raw[deep_columns].reindex(observed_pairs.sequence_right).notna().all(axis=1).to_numpy())
    observed_pairs = observed_pairs.loc[valid].reset_index(drop=True)
    check("controlled pairs have identical composition", all(sorted(a) == sorted(b) for a, b in zip(observed_pairs.sequence_left, observed_pairs.sequence_right)))
    check("controlled pairs have same source", np.array_equal(source.reindex(observed_pairs.sequence_left), source.reindex(observed_pairs.sequence_right)))
    check("controlled pairs contain no reused sequence", pd.concat([observed_pairs.sequence_left, observed_pairs.sequence_right]).is_unique)
    # A stricter dependence sensitivity: retain pairs in deterministic hash order,
    # rejecting a pair if EITHER component has appeared in an earlier kept pair.
    seen = set(); kept = []
    order = sorted(observed_pairs.index, key=lambda i: hashlib.sha256((observed_pairs.at[i, "sequence_left"] + "|" + observed_pairs.at[i, "sequence_right"]).encode()).hexdigest())
    for i in order:
        ids = {int(observed_pairs.at[i, "component_left"]), int(observed_pairs.at[i, "component_right"])}
        if seen.isdisjoint(ids):
            kept.append(i); seen.update(ids)
    pair_pops = {"all_available_pairs": observed_pairs.index,
                 "cysteine_free_length8_to39": observed_pairs.index[observed_pairs.sequence_left.str.len().between(8, 39) & ~observed_pairs.sequence_left.str.contains("C")],
                 "greedy_no_reused_component": pd.Index(kept)}
    for src in sorted(observed_pairs.source.unique()): pair_pops["source:" + src] = observed_pairs.index[observed_pairs.source.eq(src)]
    order_profiles = []; order_direction = []
    for gram, target in [("neg", "Gram-negative"), ("pos", "Gram-positive")]:
        dd = {}
        for arch in ["cnn", "lstm"]:
            col = "deep_amp_" + arch + "_" + gram + "_" + target + "_unspecified_species_and_strain"
            dd[arch] = np.log2(raw[col].reindex(observed_pairs.sequence_right).to_numpy() / raw[col].reindex(observed_pairs.sequence_left).to_numpy())
        for pop, ix in pair_pops.items():
            cnn = dd["cnn"][ix]; lstm = dd["lstm"][ix]
            for arch, change in [("CNN", cnn), ("LSTM", lstm)]:
                order_profiles.append({"target": target, "population": pop, "model": arch, "n_pairs": len(ix),
                                       "median_absolute_log2_change": np.median(abs(change)), "median_fold_change": 2 ** np.median(abs(change)),
                                       "twofold_change_n": int((abs(change) >= 1).sum()), "twofold_change_fraction": (abs(change) >= 1).mean(),
                                       "fourfold_change_fraction": (abs(change) >= 2).mean()})
            for margin in [0., .25, 1.]:
                eligible = (abs(cnn) > margin) & (abs(lstm) > margin)
                opposite = (cnn * lstm < 0) & eligible
                order_direction.append({"target": target, "population": pop, "minimum_absolute_log2_change_each": margin,
                                        "n_pairs": len(ix), "eligible_n": int(eligible.sum()), "opposite_direction_n": int(opposite.sum()),
                                        "opposite_fraction_among_eligible": fraction(opposite.sum(), eligible.sum()),
                                        "CNN_only_twofold_n": int(((abs(cnn) >= 1) & (abs(lstm) < 1)).sum()),
                                        "LSTM_only_twofold_n": int(((abs(lstm) >= 1) & (abs(cnn) < 1)).sum())})
    save(order_profiles, "DeepAMP_composition_controlled_order"); save(order_direction, "DeepAMP_order_direction_agreement")
    observed_pairs.assign(greedy_component_representative=observed_pairs.index.isin(kept)).to_csv(OUT / "DeepAMP_composition_controlled_pairs.csv", index=False)

    # Cross-check a hand-computable pairwise partition independently.
    a = mic.iloc[:, 0].le(16); b = mic.iloc[:, 1].le(16)
    d = decision_pair(mic.iloc[:, 0], mic.iloc[:, 1], 16)
    check("pairwise four cells partition cohort", sum(d[k] for k in ["both_pass_n", "left_only_n", "right_only_n", "neither_pass_n"]) == len(mic))
    check("xor agrees with marginal minus overlap identity", d["decision_disagreement_n"] == int(a.sum() + b.sum() - 2 * (a & b).sum()))
    save(CHECKS, "checks")
    summary = {"input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "n_sequences": len(raw),
               "MIC_common_n": len(mic), "classifier_complete_n": len(cp), "cutoffs_uM": CUTOFFS,
               "checks": len(CHECKS), "passed": sum(c["passed"] for c in CHECKS),
               "failed": [c for c in CHECKS if not c["passed"]], "elapsed_seconds": time.monotonic() - started,
               "pandas": pd.__version__, "numpy": np.__version__,
               "scope": "Exploratory screening decisions. Concentration cutoffs are not clinical breakpoints; classifier score thresholds are not calibrated common probabilities. Target matching follows declared export heads, not identical assays or controlled training. The DeepAMP paper reports E. coli/B. subtilis training for negative/positive heads and synthetic log-MIC labels for non-AMP examples; Gram labels do not establish general Gram-class potency."}
    (OUT / "manifest.json").write_text(json.dumps(summary, indent=2) + "\n")

    if __package__:
        from .report_figures import build_decision_figures
    else:
        from report_figures import build_decision_figures
    build_decision_figures()
    print(json.dumps(summary, indent=2)); assert all(c["passed"] for c in CHECKS)


if __name__ == "__main__":
    main()
