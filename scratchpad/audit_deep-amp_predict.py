#!/usr/bin/env python3
"""
Deep-AMP multioutput inference script for BattleAMP.

Runs all 4 model variants (LSTM/CNN x gramneg/grampos) and writes
per-variant prediction files.

Usage:
    python predict.py --input INPUT.fasta --output OUTPUT.tsv \
                      --models-dir saved_models/paper/regressor
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

MAX_LEN = 48
ALPHABET = "BCDSQKIPTFNGHLRWAVEYM-"
CHAR_TO_IDX = {c: i for i, c in enumerate(ALPHABET)}
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

VARIANTS = [
    ("LSTM_gr_neg", "predictions-lstm-gramneg.tsv"),
    ("LSTM_gr_pos", "predictions-lstm-grampos.tsv"),
    ("CNN_gr_neg",  "predictions-cnn-gramneg.tsv"),
    ("CNN_gr_pos",  "predictions-cnn-grampos.tsv"),
]

DEFAULT_VARIANT = "LSTM_gr_neg"


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_and_filter_fasta(fasta_path):
    """Read FASTA, return list of (header, sequence) tuples.

    Filters to sequences with length <= MAX_LEN and containing only
    the 20 standard amino acids.
    """
    entries = []
    header, seq_parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).upper()
                    entries.append((header, seq))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            seq = "".join(seq_parts).upper()
            entries.append((header, seq))

    valid, skipped = [], 0
    for hdr, seq in entries:
        if len(seq) < 1 or len(seq) > MAX_LEN:
            skipped += 1
            continue
        if not all(c in VALID_AA for c in seq):
            skipped += 1
            continue
        valid.append((hdr, seq))

    if skipped:
        print(f"[deep-amp] Skipped {skipped}/{len(entries)} sequences "
              f"(length or alphabet constraints)", file=sys.stderr)
    return valid


def preprocess(sequences):
    """One-hot encode sequences with padding to MAX_LEN."""
    n_chars = len(ALPHABET)
    X = np.zeros((len(sequences), MAX_LEN, n_chars), dtype=np.float32)
    for i, seq in enumerate(sequences):
        for j, aa in enumerate(seq):
            if aa in CHAR_TO_IDX:
                X[i, j, CHAR_TO_IDX[aa]] = 1.0
    return X


def write_output(sequences, predictions_log10, output_path):
    """Write predictions as TSV: sequence_id, MIC (uM), MIC_unit."""
    mic_uM = 10.0 ** predictions_log10
    rows = []
    for seq, mic in zip(sequences, mic_uM):
        rows.append({"sequence": seq, "MIC": f"{mic:.4f}", "MIC_unit": "uM"})
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


# ---------------------------------------------------------------------------
# Inference
# ---------------------------------------------------------------------------

def run_all_variants(entries, models_dir, output_base):
    """Load each model variant and produce predictions."""
    import tensorflow as tf
    tf.get_logger().setLevel("ERROR")

    headers = [h for h, _ in entries]
    sequences = [s for _, s in entries]
    X = preprocess(sequences)

    output_dir = os.path.dirname(output_base) or "."
    base_written = False

    for model_subdir, pred_filename in VARIANTS:
        model_path = os.path.join(models_dir, model_subdir)
        if not os.path.isdir(model_path):
            print(f"[deep-amp] WARNING: model not found at {model_path}, "
                  f"skipping {model_subdir}", file=sys.stderr)
            continue

        print(f"[deep-amp] Loading {model_subdir}...", file=sys.stderr)
        model = tf.keras.models.load_model(model_path, compile=False)
        preds = model.predict(X, verbose=0).flatten()

        variant_output = os.path.join(output_dir, pred_filename)
        write_output(sequences, preds, variant_output)
        print(f"[deep-amp] Wrote {len(preds)} predictions to {variant_output}",
              file=sys.stderr)

        # Write base output (default variant = LSTM gram-neg)
        if model_subdir == DEFAULT_VARIANT and not base_written:
            write_output(sequences, preds, output_base)
            base_written = True

    if not base_written:
        print("[deep-amp] WARNING: default variant not found, "
              "base output not written", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Deep-AMP multioutput inference")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Base output TSV path")
    parser.add_argument("--models-dir", required=True,
                        help="Directory containing the 4 SavedModel folders")
    args = parser.parse_args()

    entries = read_and_filter_fasta(args.input)
    if not entries:
        print("[deep-amp] No valid sequences after filtering.", file=sys.stderr)
        sys.exit(1)

    print(f"[deep-amp] {len(entries)} sequences passed filters", file=sys.stderr)
    run_all_variants(entries, args.models_dir, args.output)


if __name__ == "__main__":
    main()
