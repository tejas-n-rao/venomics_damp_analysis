# AMPredictor

Fork of [ruihan-dong/AMPredictor](https://github.com/ruihan-dong/AMPredictor), adapted
for integration with the
[battleamp-snakemake](https://github.com/szczurek-lab/battleamp-snakemake) benchmarking
pipeline.

## Supported tasks

AMP regression (predicts minimum inhibitory concentration in uM).

## Reference

Dong, R., Liu, R., Liu, Z. et al. Exploring the repository of de novo designed
bifunctional antimicrobial peptides through deep learning. *eLife* (2025).
https://doi.org/10.7554/eLife.97330.3

## Requirements

- Python 3.11
- conda (for environment creation by the pipeline)
- NVIDIA GPU (ESM-1b and GCN inference)

The model architecture (ESM-1b embeddings fed into a graph convolutional network) and
pretrained weights are unchanged. The following issues were fixed and files added:

**Inference optimization:** The original pipeline wrote ESM embeddings, contact maps,
molecular graphs, and fingerprints to disk as intermediate files, then read them back.
This was replaced with an in-memory pipeline: proper ESM batching, chunked processing,
no intermediate disk I/O, and precomputed fingerprints. This significantly reduces
inference time on large datasets. The original model weights, training data, and GCN architecture code are preserved as-is.


## Installation

```bash
conda create -n ampredictor python=3.11.9
conda activate ampredictor
sh setup.sh
```

Test whether everything works:

```bash
sh inference.sh sample.fasta results.tsv
```

## Usage within the pipeline

This repository is included as a git submodule in battleamp-snakemake:

```bash
cd battleamp-snakemake
git submodule add git@github.com:szczurek-lab/BattleAMP-AMPredictor.git models/ampredictor
```

The pipeline handles environment creation, inference, and evaluation automatically.

## Notes

- Maximum sequence length is 65 amino acids. Longer sequences are skipped by the
  pipeline's prefilter step (coverage is reported alongside metrics).
- The model was trained exclusively on antimicrobial peptides with known MIC values. It
  has no concept of "inactive" and always predicts a MIC in a plausible range. Binary
  classification via thresholding is handled by the benchmark pipeline's
  regressor-to-classification conversion, but ranking metrics (AUROC, Spearman) are more
  informative than binary metrics for this model.
- ESM-1b weights (~2.5 GB) are downloaded once during setup and cached locally.
- Output MIC values are in micromolar (uM).

## License

Same as the original AMPredictor repository.