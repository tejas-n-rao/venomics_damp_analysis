"""Record the model-scope review and a selected-output sensitivity for DeepAMP.

python analysis_review/review_mic_model_scope.py
Paper annotations supplement rather than overwrite upstream export metadata.
Training-data provenance of individual checkpoints is not independently rebuilt.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_review/outputs/model_scope_review"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    path = ROOT / "df_damp_195625x57.tsv"
    columns = ["deep_amp_" + arch + "_" + gram + "_" + target + "_unspecified_species_and_strain"
               for arch in ["cnn", "lstm"] for gram, target in [("neg", "Gram-negative"), ("pos", "Gram-positive")]]
    raw = pd.read_csv(path, sep="\t", usecols=["sequence"] + columns).set_index("sequence", verify_integrity=True)
    pairs = pd.read_csv(ROOT / "analysis_review/outputs/screening_decision_trends/DeepAMP_composition_controlled_pairs.csv")
    rows = []
    for gram, target in [("neg", "Gram-negative"), ("pos", "Gram-positive")]:
        left = []; right = []
        for arch in ["cnn", "lstm"]:
            col = "deep_amp_" + arch + "_" + gram + "_" + target + "_unspecified_species_and_strain"
            left.append(raw[col].reindex(pairs.sequence_left).to_numpy())
            right.append(raw[col].reindex(pairs.sequence_right).to_numpy())
        ll = np.array(left).T; rr = np.array(right).T
        changes = np.abs(np.log2(rr / ll))
        for cap in [32., 128., 10 ** 3.5]:
            eligible = (ll <= cap).all(axis=1) & (rr <= cap).all(axis=1)
            for j, arch in enumerate(["CNN", "LSTM"]):
                rows.append({"export_target_label": target, "paper_training_species": "E. coli" if gram == "neg" else "B. subtilis",
                             "maximum_prediction_uM_each_sequence_and_model": cap, "model": arch,
                             "original_pair_n": len(pairs), "selected_pair_n": int(eligible.sum()),
                             "twofold_change_n": int((changes[eligible, j] >= 1).sum()),
                             "twofold_change_fraction": (changes[eligible, j] >= 1).mean(),
                             "interpretation": "Conditions on all four predictions within this target; descriptive, not independent biological validation or a causal test of pseudo-label effects."})
    pd.DataFrame(rows).to_csv(OUT / "DeepAMP_order_low_prediction_sensitivity.csv", index=False)
    meta = pd.DataFrame(json.loads((ROOT / "scratchpad/damp_metadata.json").read_text()))
    sources = {
        "DeepAMP_paper": "https://www.nature.com/articles/s41467-023-42434-9",
        "MBC_paper": "https://journals.asm.org/doi/10.1128/msystems.00345-23",
        "AMPredictor_paper": "https://elifesciences.org/articles/97330",
        "APEX_1_1_paper": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12408343/",
    }
    review = {"review_date": "2026-11-09", "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "metadata_sha256": hashlib.sha256((ROOT / "scratchpad/damp_metadata.json").read_bytes()).hexdigest(),
              "source_urls": sources,
              "DeepAMP_paper_reported_training": {"Gram-negative_export": "4,619 unique AMP sequences with E. coli MIC measurements",
                                                  "Gram-positive_export": "4,089 unique AMP sequences with B. subtilis MIC measurements",
                                                  "non_AMP_label_log10_MIC": 4,
                                                  "checkpoint_specific_training_provenance_rebuilt": False},
              "MBC_training_target": "E. coli measurements pooled across conditions; paper reports averaging repeated measurements. Exact averaging stage not independently reconstructed.",
              "AMPredictor_training_scope": "Paper reports 3,280 E. coli AMP sequences, no cysteine and lengths below 40; accepted parser lengths are not the training domain.",
              "APEX_1_1_inactive_pseudolabel_uM": 512,
              "pinned_export_revisions": meta[(meta.model.isin(["apex_1_0", "apex_1_1", "mbc_attention", "ampredictor"]) | meta.model.str.startswith("deep_amp"))][["model", "revision"]].drop_duplicates().to_dict("records"),
              "policy": "Keep native column labels and numerical tables unchanged; correct prose and display annotations. No cross-species pooling, biological error ranking, or causal architecture attribution.",
              "sensitivity_note": "Output-conditioned low-MIC subsets cannot establish that synthetic inactive labels cause or do not cause model differences."}
    (OUT / "review_record.json").write_text(json.dumps(review, indent=2) + "\n")
    print(pd.DataFrame(rows)[["export_target_label", "maximum_prediction_uM_each_sequence_and_model", "model", "selected_pair_n", "twofold_change_fraction"]].to_string(index=False))


if __name__ == "__main__":
    main()
