# DAMP analysis

The main findings are detailed in [final report](reports/FINAL_ANALYSIS.md). You could also use the accompanying html in your local browser. However, if you'd want to know more for a certain subsection, the markdown links to "companion reports" (find here: [reports/](reports/)). These expound upon the specifc analysis in the final report. 

## Focused reports

<details markdown="1">
<summary>Show focused reports and supporting analyses</summary>

| Final-report section | Detailed report |
|---|---|
| 1. Different models choose different peptides | [Shortlists and robustness](reports/01_shortlists.md) |
| 2. A concentration cutoff changes the screening decision | [Concentration rules and MIC context](reports/02_concentration_rules.md) |
| 3. Differences remain when the declared target is matched | [Matched targets, model training, and conversions](reports/03_matched_targets.md) |
| 4. General activity classifiers often split their decisions | [Classifier groups and sequence analysis](reports/04_classifier_sequences.md) |
| 5. Sequence order affects models differently even for the same target | [Sequence responses](reports/05_sequence_responses.md) |
| 6. Venomics gets more unanimous positive votes, while MIC shortlists differ | [Venomics source comparison](reports/06_venomics.md) |
| 7. What the results support | [Interpretation and experimental design](reports/07_interpretation.md) |

[Source comparisons](reports/SOURCE_COMPARISONS.md) follows the main questions separately in all eight collections, with checks for related variants and different denominators. Venomics is reported first; other sources stand out on different dimensions.

[Preliminary analyses](reports/PRELIMINARY_ANALYSES.md) collects agreement-definition profiles, wider output patterns, and provisional experimental-data recovery. [Reproducibility](reports/REPRODUCIBILITY.md) holds calculation and inference checks. [Methods and the analysis map](reports/METHODS.md) retain definitions, detailed controls, and the notebook-to-report guide.

</details>

## Reproducibility

<details markdown="1">
<summary>Show notebooks, reproducibility checks, and rebuild instructions</summary>

The [reviewed notebook](damp_secondary_analysis.ipynb) and [early EDA](damp_exploratory_analysis.ipynb) retain their code and saved outputs. Scientific tables remain in [analysis_outputs](analysis_outputs/) and [analysis_review/outputs](analysis_review/outputs/). Pinned source/model snapshots remain with their provenance inputs.


### Rebuild the figures and reports

Use the existing Python environment (`python`, rather than the older `python3` in this workspace):

```bash
python analysis_review/plot_shortlist_upset.py
python analysis_review/report_figures.py
python analysis_review/build_readable_report.py --all
```

The shared plot module reads saved analysis tables. Rebuilding reports does not fit models or change scientific results. A single report can also be rendered with `--source reports/04_classifier_sequences.md`.

### Reproduce numerical analyses

```bash
python analysis_review/review_analysis.py
python analysis_review/reproduce_notebook.py
python analysis_review/curate_archaeasin_labels.py
python analysis_review/stress_test_findings.py
python analysis_review/explore_additional_trends.py
python analysis_review/sequence_disagreement.py
python analysis_review/similar_mic_different_ranks.py
python analysis_review/screening_decision_trends.py
python analysis_review/review_mic_model_scope.py
python analysis_review/classifier_sequence_analysis.py
python analysis_review/classifier_rule_sensitivity.py
python analysis_review/venomics_subgroup.py
python analysis_review/verify_venomics_subgroup.py
python analysis_review/source_stratified_analysis.py
python analysis_review/verify_source_stratified_analysis.py
```

The notebook rerun needs the cached DeepAMP environment and weights. Experimental curation uses the downloads listed in its [manifest](analysis_review/outputs/archaeasin_curation_manifest.json). Later sequence scripts consume the stress-test component map and additional-trends composition pairs; the classifier extension also consumes the exact-substitution catalogue. [Support-file notes](scratchpad/README.md) describe retained inputs and legacy references.

Duplicate summaries, accessible versions, and superseded report snapshots were removed after their unique material was merged into the canonical reports. The [consolidation manifest](analysis_review/markdown_consolidation_manifest.json) records routes, removed filenames, and checksums. The earlier [code-cleanup manifest](analysis_review/cleanup_manifest.json) records prior support-file removals. Scientific inputs, numerical tables, notebook code/outputs, and pinned source/model snapshots are preserved.

</details>
