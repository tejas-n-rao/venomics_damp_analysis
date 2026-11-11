# How the calculations and inference routes were checked

Reproducibility establishes what was calculated from the saved data. Biological validation asks whether predictions match the intended experiments. This report collects implementation evidence so the [final analysis](FINAL_ANALYSIS.md) can concentrate on scientific findings.

## 1. Original calculations and independent recomputation

A **calculation rerun** executes the analysis again from saved predictions; it does not regenerate those predictions from model weights. The comparison checks table names and columns as well as numerical values.

An isolated rerun completed the original reviewed notebook's **20 code cells** and reproduced **72 generated CSV tables**. Numeric comparisons allowed a difference no larger than 10⁻⁸ plus 10⁻⁶ times the reference value (an absolute floor plus a size-dependent tolerance); schema and text matched exactly. Later scripts have their own checks and are not presented as part of that original executed notebook.

An independent raw-data review used a separate implementation to check the calculations. A **hash** is a checksum identifying a file's contents; matching hashes help detect unintended changes. The review checked input/source hashes, sequence alignment, finite positive MIC outputs, group membership, matching rules, shortlists, and arithmetic identities. Its **99 checks** passed. Repeating a calculation or asserting an identity supports implementation consistency; it is not independent biological replication. [Notebook reproduction](../analysis_review/outputs/notebook_reproduction.json), [independent checks](../analysis_review/outputs/independent_checks.csv), and [review summary](../analysis_review/outputs/review_summary.json).

## 2. A sampled DeepAMP inference route matched the saved table

**Inference replay** runs saved model weights on sequences and compares the resulting predictions with the exported table. **Padding** fills unused positions so variable-length peptides can enter a fixed-size model input; different filling rules can change the predictions. A log2 replay/export ratio of zero means an exact match, +1 means replay is twice the export, and −1 means half. The absolute value measures the size of the discrepancy in either direction.

Four saved models were evaluated on **40 selected sequences**: 32 seeded random eligible entries and eight extremes. Across 160 predictions, the token-padding route reproduced the TSV within **6.2 × 10⁻⁶ absolute log2 units**. This is a numerical reproduction check on the selected sequences.

An alternative zero-padding route with the same pinned weights changed predictions by up to **8.52 log2 units**, roughly **367-fold**, in this selected sample. This establishes the importance of preprocessing and identifies a sampled route consistent with the exports. It does not show that every saved value used that route or reconstruct every model's inference. Full APEX, MBC, and AMPredictor weight replays remain outstanding. [Replay summary](../analysis_outputs/deepamp_replay_summary.csv) and [provenance](../analysis_outputs/replay_manifest.json).

## 3. Later extensions have explicit definitions and checks

An **implementation check** is a software assertion about a defined calculation, such as whether exclusive overlap counts sum to the union or whether all selected pairs have the stated sequence relationship. Several checks can use the same data.

The stress tests record 32 checks, additional trends 25, sequence-disagreement extension 15, similar-MIC subgroup 12, screening-decision script 126, and exact UpSet calculation eight. These counts describe assertions across definitions and thresholds, not quantities of independent scientific evidence. Their scripts, dependencies, and manifests are linked in the [analysis map](METHODS.md#analysis-map).

The shortlist extension also checks the top-2% and top-5% intersections: every list size, partition of the union, reconstruction of total selections, independent set intersections, and untied boundary. Spearman correlation measures consistency of peptide ordering. Its estimates are checked against correlation of average tied ranks (equal predictions receive the average of their occupied ranks), and the full-cohort values reproduce the earlier review. [Budget and correlation checks](../analysis_review/outputs/shortlist_overlap/budget_extension_manifest.json).

The 16 µM cutoff UpSet separately checks all eight membership combinations, including joint failures. Seven assertions confirm the cohort, partitions, each model's total, independent set calculations, and agreement with the earlier cutoff analysis. All pass. [Cutoff-overlap manifest](../analysis_review/outputs/shortlist_overlap/cutoff_16uM/manifest.json).

The classifier-sequence extension checks its 195,620-sequence denominator, vote totals against the separately saved 16-pattern table, group counts, residue accounting, zero length differences after exact-length matching, and validity of the 41,253 one-substitution pairs. It preserves unmapped sequences explicitly, rather than inventing group assignments. [Classifier code](../analysis_review/classifier_sequence_analysis.py), [checks](../analysis_review/outputs/classifier_sequences/checks.csv), and [manifest](../analysis_review/outputs/classifier_sequences/manifest.json).

Here, K/R means the combined lysine/arginine fraction of a sequence. A neighbor pair differs at one amino-acid position; a related group joins observed substitution and terminal-extension links. A **source/length-standardized estimate** compares sequences only within the same source collection and exact length before averaging.

An independent raw-data calculation also recovers the 32,423 same-group neighbor pairs, each pair's changed classifier decisions, and the source/length-standardized K/R estimate. Explicit comparator sets for 24 spread-out pairs verify the ordinary and outside-component reference formulas. It confirms 73 unanimous-positive/unanimous-negative switches across 11 groups, with 62 in the largest three. These checks use the same saved predictions and are not independent outcome evidence. [Independent verification](../analysis_review/outputs/classifier_sequences/independent_verification.json).

The two-rule classifier extension uses a separate [script](../analysis_review/classifier_rule_sensitivity.py) to repeat group counts, source/length-matched sequence contrasts, and neighbor/reference agreement at 0.3, 0.5, and 0.7. Its 22 checks include agreement with the earlier saved profiles, monotonic vote changes, unchanged binary AMP-CLIP decisions, sequence-derived K/R and D/E accounting, and explicit reference sets for 24 anchors per rule. The underlying pairs and score-complete population stay fixed. [Results and input hashes](../analysis_review/outputs/classifier_sequences/rule_sensitivity/manifest.json).

The Venomics extension preserves the primary 194,947-sequence cohort and all 4,618 Venomics rows. Seven primary checks verify cohort coverage, classifier partitions, the 401 related groups, cutoff counts, exact C-presence matching, and bounded rates. A [separate raw-data script](../analysis_review/verify_venomics_subgroup.py) passes 50 additional checks, including feature counts, source/group-weighted rates, every 5% and cutoff overlap cell, explicit matched-category arithmetic, and all 24 source correlations against the prior review. These are arithmetic checks on the same predictions, not independent biological validation. [Analysis manifest](../analysis_review/outputs/venomics/manifest.json) and [verification results](../analysis_review/outputs/venomics/independent_verification.json).

The all-source extension adds 10 primary checks and 226 [separate arithmetic checks](../analysis_review/outputs/source_comparisons/independent_verification.json). These cover within-source shortlists, a prediction-independent representative per group, matched-target decision cells, length-matched amino-acid contrasts, and source-specific order responses. The source/pair counts and different eligible populations are explicit in the [manifest](../analysis_review/outputs/source_comparisons/manifest.json). Repeating checks on the same predictions verifies the calculation, not biological generalization.

## 4. Reporting changes preserve scientific inputs

The section-specific reports embed figures and link their underlying data. The report renderer rebuilds HTML/PDF from Markdown; it does not fit models or rerun scientific analyses. The shared plotting module reads the saved result tables. The original prediction TSV and notebook code/saved outputs are preserved. Numerical tables from earlier analyses are checked for unintended changes when adding a new result branch.

The [scientific audit](METHODS.md) retains methodological history, claim corrections, and full metric definitions. The [validation findings](07_interpretation.md#2-existing-assay-matches-do-not-identify-a-winner) and [section 7](07_interpretation.md) explain why this computational evidence does not yet establish an accuracy winner.
