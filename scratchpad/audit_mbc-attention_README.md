# MBC-Attention

Fork of [jieluyan/MBC-Attention](https://github.com/jieluyan/MBC-Attention), adapted for
integration with the
[battleamp-snakemake](https://github.com/szczurek-lab/battleamp-snakemake) benchmarking
pipeline.

## Supported tasks

AMP regression (predicts minimum inhibitory concentration for E. coli).

## Reference

Yan J, Zhang B, Zhou M, Campbell-Valois F, Siu SWI.2023.A deep learning method for predicting the minimum inhibitory concentration of antimicrobial peptides against Escherichia coli using Multi-Branch-CNN and Attention. mSystems8:e00345-23.https://doi.org/10.1128/msystems.00345-23

## Model overview

MBC-Attention is a multi-branch CNN with attention that predicts minimum inhibitory
concentration (MIC) values for antimicrobial peptides against E. coli. Input sequences
are featurized using 14 sequence-derived descriptors from iFeature (PseKRAAC and
QSOrder), reshaped into 2D feature images, and processed by parallel convolutional
branches with an attention mechanism.

The model predicts MIC in micromolar (uM) using log10-transformed values. Lower
predicted values indicate higher predicted antimicrobial activity.

Sequence length constraints: 5 to 60 amino acids (imposed by the QSOrder descriptor).


## Requirements

- conda (for environment creation by the pipeline)
- GPU optional (TF auto-detects; CPU fallback supported)
- Model weights: `model/whole_train.mdl/` (TF SavedModel directory, included)
- Feature data: `data/EC.csv` (included)

## Installation

```bash
conda env create -f environment.yaml
conda activate mbc-attention
sh setup.sh
```

## Usage within the pipeline

This repository is included as a git submodule in battleamp-snakemake:

```bash
cd battleamp-snakemake
git submodule add git@github.com:szczurek-lab/BattleAMP-mbc-attention.git models/mbc-attention
```

## Standalone usage

```bash
sh inference.sh input.fasta predictions.tsv
```

Output columns: `sequence`, `Prediction`, `Probability_score`. `Prediction` and
`Probability_score` are the predicted MIC in micromolar. Lower values indicate higher
predicted antimicrobial potency.

## Notes

- Sequences outside the 5-60 AA length range are silently excluded by the model's
  internal preprocessing (config.py min/max_seq_length).
- The model uses iFeature for descriptor computation, which is bundled in the
  `ifeature/` directory.
- PyTorch is required as a transitive import (`tools/__init__.py`) but is not used in
  the inference path. Only TF/Keras is used for prediction.
- The LD_LIBRARY_PATH fix in inference.sh handles TF's pip-installed NVIDIA libraries
  for GPU detection.

## License

Same as the original MBC-Attention repository.