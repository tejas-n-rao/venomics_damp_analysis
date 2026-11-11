"""Recover a transparently provisional, chemistry-audited archaeasin diagnostic.

Inputs are the published Fig 2 / Extended Data Fig 4 workbooks and the Mendeley
v3 synthesis certificates. The experimental workbook has peptide numbers, not
sequences. We therefore infer a candidate sequence using its 11-target published
prediction fingerprint, then require an exact synthesis-certificate sequence.
Measured MICs are NEVER used to select the sequence match. Number-to-sequence
linkage remains inferred and should be author-confirmed before a benchmark claim.
"""
from pathlib import Path
import hashlib
import io
import json
import os
import re
import tempfile
import zipfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "venomics-mpl"))
os.environ["MPLBACKEND"] = "Agg"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "analysis_review/inputs"
OUT = ROOT / "analysis_review/outputs"
BASE = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/"
URLS = {
    "archaeasin_figure2_source.xlsx": BASE + "41564_2025_2061_MOESM8_ESM.xlsx",
    "archaeasin_extended4_source.xlsx": BASE + "41564_2025_2061_MOESM14_ESM.xlsx",
    "archaeasin_sequences.csv": BASE + "41564_2025_2061_MOESM4_ESM.csv",
    "archaeasin_supplement.pdf": BASE + "41564_2025_2061_MOESM1_ESM.pdf",
    "archaeasin_certificates.zip": "https://data.mendeley.com/public-files/datasets/d8yzgtdrcp/files/b1307377-6e4d-4f63-b2c1-da070f075949/file_downloaded",
}


def main():
    OUT.mkdir(exist_ok=True)
    archive_hash = hashlib.sha256((IN / "archaeasin_certificates.zip").read_bytes()).hexdigest()
    assert archive_hash == "982cd9a9f6d497c14c67a1b72ea25106cc865743ae63f2e1275d73bad37ef0b3"
    certificate_rows = []
    with zipfile.ZipFile(IN / "archaeasin_certificates.zip") as archive:
        for name in archive.namelist():
            if not name.lower().endswith(".pdf") or name.startswith("__MACOSX"):
                continue
            text = PdfReader(io.BytesIO(archive.read(name))).pages[0].extract_text()
            match = re.search(r"Peptide Sequence\s*:\s*([ACDEFGHIKLMNPQRSTVWY]+)", text)
            assert match, name
            certificate_rows.append({"certificate_file": name, "sequence": match.group(1),
                                     "salt_form": "TFA" if "Salt Form:TFA" in text else "unresolved",
                                     "terminal_modifications": "not explicitly specified in extracted certificate text"})
    certificates = pd.DataFrame(certificate_rows).set_index("sequence", verify_integrity=True)
    assert len(certificates) == 80
    certificates.to_csv(OUT / "archaeasin_certificate_sequences.csv")

    data = pd.read_csv(ROOT / "df_damp_195625x57.tsv", sep="\t").set_index("sequence", verify_integrity=True)
    meta = pd.DataFrame(json.loads((ROOT / "scratchpad/damp_metadata.json").read_text())).set_index("column")
    targets = list(meta.index[meta.model.eq("apex_1_1")])
    workbook = IN / "archaeasin_extended4_source.xlsx"
    sheets = pd.ExcelFile(workbook).sheet_names[1:12]
    assert len(targets) == len(sheets) == 11
    exported = np.log2(data[targets].to_numpy())
    # Native output indices and explicit worksheet strain names verify the order.
    expected_tokens = ["ATCC19606", "ATCC11775", "AIC221", "AIC222", "ATCC13883",
                       "PAO1", "PA14", "ATCC12600", "ATCCBAA1556", "ATCC700802", "ATCC700221"]
    assert all(token in sheet for token, sheet in zip(expected_tokens, sheets))
    assert meta.loc[targets, "native_output_index"].tolist() == list(range(11))
    predictions = np.column_stack([pd.read_excel(workbook, sheet).iloc[:, 0].to_numpy(dtype=float) for sheet in sheets])
    source_labels = np.column_stack([pd.read_excel(workbook, sheet).iloc[:, 1].to_numpy(dtype=float) for sheet in sheets])
    measurements = pd.read_excel(IN / "archaeasin_figure2_source.xlsx", sheet_name="Fig 2a", index_col=0)
    assert measurements.shape == (80, 11)
    assert list(measurements.index) == [f"Archaeasin-{i}" for i in range(1, 81)]
    numeric = measurements.notna().to_numpy()
    assert np.allclose(np.log2(measurements.to_numpy(dtype=float)[numeric]), source_labels[numeric])
    assert np.all(source_labels[~numeric] == 9)  # 512 is the publication's placeholder.
    shortlist = set(pd.read_csv(IN / "archaeasin_sequences.csv").iloc[:, 0])
    mappings = []
    for row_index, vector in enumerate(predictions):
        max_error = np.max(abs(exported - vector), axis=1)
        nearest = np.argsort(max_error)[:2]
        sequence = data.index[nearest[0]]
        # Tight multi-target match and a separated runner-up, required for every row.
        assert max_error[nearest[0]] < .01
        assert max_error[nearest[1]] > .05
        leave_target_out_agrees = all(
            np.argmin(np.max(abs(np.delete(exported, j, axis=1) - np.delete(vector, j)), axis=1)) == nearest[0]
            for j in range(11)
        )
        assert leave_target_out_agrees
        mappings.append({"peptide_id": f"Archaeasin-{row_index + 1}", "sequence": sequence,
                         "linkage_method": "11-target prediction fingerprint; experimental labels unused",
                         "max_log2_prediction_difference": max_error[nearest[0]],
                         "runner_up_max_log2_difference": max_error[nearest[1]],
                         "all_11_leave_target_out_matches_agree": leave_target_out_agrees,
                         "in_published_265_sequence_shortlist": sequence in shortlist,
                         "exact_sequence_in_synthesis_certificates": sequence in certificates.index})
    mapping = pd.DataFrame(mappings)
    assert mapping.sequence.is_unique
    mapping.to_csv(OUT / "archaeasin_provisional_id_map.csv", index=False)
    # Exclude mismatched certified chemistry instead of silently accepting a terminal variant.
    missing_certified = sorted(set(certificates.index) - set(data.index))
    assert missing_certified == ["VKWKKLIKGFVINKF"]
    assert mapping.loc[~mapping.exact_sequence_in_synthesis_certificates, "sequence"].tolist() == ["KWKKLIKGFVINKF"]
    accepted = mapping.set_index("peptide_id")
    accepted = accepted[accepted.exact_sequence_in_synthesis_certificates].copy()
    accepted["MIC_uM"] = measurements.iloc[:, 1].reindex(accepted.index)
    accepted["published_evaluation_placeholder_uM"] = np.exp2(source_labels[:, 1])[mapping.exact_sequence_in_synthesis_certificates.to_numpy()]
    accepted["label_status"] = np.where(accepted.MIC_uM.notna(), "numeric", "blank_Fig2; coded_512_in_ED4")
    accepted["target"] = "E. coli ATCC 11775"
    accepted["chemistry_limit"] = "exact certificate base sequence and TFA salt; terminal chemistry independently unresolved"
    db = pd.read_csv(ROOT / "scratchpad/replay/APEXDB.csv")
    amp = set()
    for line in (ROOT / "scratchpad/replay/ampredictor_train.txt").read_text().splitlines():
        fields = line.replace("\\t", "\t").split("\t")
        if len(fields) >= 2:
            amp.add(fields[1].strip())
    accepted["exact_APEXDB_overlap"] = accepted.sequence.isin(db.Sequence)
    accepted["exact_AMPredictor_training_overlap"] = accepted.sequence.isin(amp)
    models = {"APEX 1.0": "apex_1_0_E._coli_ATCC11775", "APEX 1.1": "apex_1_1_E._coli_ATCC_11775",
              "MBC-Attention": "mbc_attention_E._coli_unspecified_strain", "AMPredictor": "ampredictor_E._coli_unspecified_strain"}
    for model, column in models.items():
        accepted[model] = data.loc[accepted.sequence, column].to_numpy()
    accepted.to_csv(OUT / "archaeasin_provisional_ecoli_labels.csv")
    metrics = []
    for model in models:
        block = accepted.dropna(subset=["MIC_uM", model])
        error = np.log2(block[model] / block.MIC_uM)
        metrics.append({"model": model, "n_numeric": len(block), "MAE_log2_numeric_only": error.abs().mean(),
                        "median_signed_error_log2": error.median(), "within_twofold_numeric_only": error.abs().le(1).mean(),
                        "spearman_numeric_only": block[model].corr(block.MIC_uM, method="spearman"),
                        "n_source_coded_inactive": int(accepted.MIC_uM.isna().sum()),
                        "inactive_predictions_le_64_fraction": accepted.loc[accepted.MIC_uM.isna(), model].le(64).mean(),
                        "inactive_interpretation": "conditional diagnostic: source blanks interpreted as >64; never exact MIC=512"})
    pd.DataFrame(metrics).to_csv(OUT / "archaeasin_provisional_diagnostics.csv", index=False)
    # Demonstrate the impact of publication label coding using source values alone;
    # no sequence linkage is needed for this calculation.
    coding_rows = []
    for j, target in enumerate(measurements.columns):
        numeric_j = numeric[:, j]
        for subset, mask in [("all_including_512_placeholders", np.ones(80, dtype=bool)), ("numeric_observations_only", numeric_j)]:
            a = pd.Series(predictions[mask, j]); b = pd.Series(source_labels[mask, j])
            coding_rows.append({"target": target, "subset": subset, "n": int(mask.sum()),
                                "pearson_log2": a.corr(b), "spearman": a.corr(b, method="spearman")})
    pd.DataFrame(coding_rows).to_csv(OUT / "archaeasin_publication_label_coding.csv", index=False)
    manifest = {"article": "https://doi.org/10.1038/s41564-025-02061-0", "retrieved_date_India": "2026-11-06",
                "files": [{"name": name, "url": url, "sha256": hashlib.sha256((IN / name).read_bytes()).hexdigest()} for name, url in URLS.items()],
                "n_source_peptides": 80, "n_certificate_base_sequence_matches": len(accepted),
                "n_numeric_Ecoli_after_chemistry_exclusion": int(accepted.MIC_uM.notna().sum()),
                "n_blank_source_coded_inactive": int(accepted.MIC_uM.isna().sum()),
                "n_exact_APEXDB_overlap": int(accepted.exact_APEXDB_overlap.sum()),
                "n_exact_AMPredictor_training_overlap": int(accepted.exact_AMPredictor_training_overlap.sum()),
                "mismatched_prediction_sequence": "KWKKLIKGFVINKF", "unmatched_certified_sequence": "VKWKKLIKGFVINKF",
                "excluded_peptide": "Archaeasin-31", "excluded_numeric_Ecoli_MIC_uM": 32,
                "id_linkage_status": "inferred from unique prediction fingerprints, not author-confirmed identifiers",
                "evaluation_status": "selected-cohort provisional diagnostic; not an independent multi-model benchmark"}
    (OUT / "archaeasin_curation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k:v for k,v in manifest.items() if k != "files"}, indent=2))
    print(pd.DataFrame(metrics).to_string(index=False))


if __name__ == "__main__":
    main()
