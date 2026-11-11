# Analysis support files

The executed [reviewed notebook](../damp_secondary_analysis.ipynb) and numerical scripts remain the reproducible sources. The current [report index](../README.md) gives the reading and rebuild paths.

## Inputs and provenance

- `damp_metadata.json` is the metadata used by the notebook. `upstream_metadata.json` is its byte-identical downloaded provenance snapshot; both names remain because the analysis and source references distinguish their roles.
- `source_sequences.tsv` maps exact sequences to source labels, verified against the pinned upstream source tables. `marlys_source.tsv` is the materialized append table.
- `DAMP/` contains the cached upstream checkout, excluded from Git. Its revision is recorded in the analysis manifest.
- `audit_*` files preserve model README/code snapshots at recorded revisions. They support the [model-comparison report](../reports/03_matched_targets.md#model-targets-training-and-conversions) and are evidence rather than unused application code.
- `replay/` contains downloaded development-corpus lists, replay metadata, the isolated environment, and model weights. These support sampled inference and overlap checks.

The input TSV is excluded from Git but required for numerical reproduction. Its hash is recorded in the analysis manifests. Cached source files, experimental inputs, and numeric results are retained.

## Runnable diagnostics

- `replay_deepamp.py` uses `replay/venv/bin/python` and cached SavedModels. Its weight/source records are retained in `replay/downloads.json` and the replay manifest.
- `mic_variability.py` computes tolerance and paired-label summaries; notebook section 21 executes it.
- `experimental_validation.py` and `training_overlap.py` provide standalone diagnostic calculations. Their exposed assay sample is not an independent accuracy benchmark.
- `sequence_similarity.py`, `low_agreement.py`, and `high_low_disagreement.py` retain the earlier numerical investigations. The reviewed scripts in `analysis_review/` add target, domain, and dependence checks rather than silently replacing their saved results.

The original Colab probes (`1_DATA_PROBE.py`, `3_apex_verification.py`, and `3b_MARLYS_COMPARISONS.py`) are retained as unique historical reference material. They contain notebook-specific syntax and are not part of the local rerun workflow.

Old scripts that built or appended notebook/report prose were retired: rerunning them would overwrite the reviewed wording or duplicate sections. Unreferenced GitHub directory listings and temporary run logs were also removed. The [cleanup manifest](../analysis_review/cleanup_manifest.json) records the removals; the current reports are rendered with one shared builder.
