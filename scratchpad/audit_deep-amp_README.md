# Deep-AMP

Adapter for [Deep-AMP](https://github.com/amirpandi/Deep_AMP), integrated with the
[battleamp-snakemake](https://github.com/szczurek-lab/battleamp-snakemake) benchmarking
pipeline.

## Supported tasks

AMP regression (predicts MIC via classification-through-regression).

## Reference

Pandi, A., Adam, D., Zare, A. et al. Cell-free biosynthesis combined with deep
learning accelerates de novo-development of antimicrobial peptides. *Nature
Communications* 14, 7197 (2023). https://doi.org/10.1038/s41467-023-42434-9

## Model variants

| Variant | Architecture | Training data |
|---|---|---|
| `deep-amp-lstm-gramneg` | LSTM(500) | Gram-negative MIC |
| `deep-amp-lstm-grampos` | LSTM(500) | Gram-positive MIC |
| `deep-amp-cnn-gramneg` | Conv1D | Gram-negative MIC |
| `deep-amp-cnn-grampos` | Conv1D | Gram-positive MIC |

All four run in a single inference call (multi-output pattern).


## Requirements

- Python 3.10
- conda (for environment creation by the pipeline)
- GPU not required (TensorFlow SavedModel, runs on CPU)
The model architecture and pretrained weights are unchanged. Output MIC values are in
micromolar (via `10^predicted_log10_MIC`).

## Installation

```bash
conda create -n deep-amp python=3.10
conda activate deep-amp
sh setup.sh
```

## Constraints

- Maximum sequence length: 48 amino acids.
- Only the 20 standard amino acids are accepted.

## Notes

- The Deep-AMP training data is also available as a BattleAMP benchmark dataset. Note
  that Deep-AMP models were trained on this data, so results on this dataset reflect
  memorization rather than generalization. This is documented for transparency.

## License

Same as the original Deep-AMP repository.