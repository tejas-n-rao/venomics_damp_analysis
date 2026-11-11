import json
import os
import sys
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Batch
from torch_geometric import data as DATA
from torch.utils.data.dataloader import default_collate

import esm
from Bio import SeqIO
from rdkit.Chem import AllChem as Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')

from AMPredictor import GNNPredictor

MAX_LENGTH = 65
ESM_LAYER = 33
ESM_DIM = 1280
CONTACT_THRESHOLD = 0.5
FP_RADIUS = 3
FP_NBITS = 2048


# ---------------------------------------------------------------------------
# 1. Read and filter FASTA
# ---------------------------------------------------------------------------

def read_fasta(fasta_path):
    """Read FASTA, return (ids, sequences) for sequences <= MAX_LENGTH."""
    ids, seqs = [], []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        if len(seq) <= MAX_LENGTH:
            ids.append(str(record.id))
            seqs.append(seq)
    return ids, seqs


# ---------------------------------------------------------------------------
# 2. Fingerprint features (precomputed lookup)
# ---------------------------------------------------------------------------

def precompute_aa_fingerprints(smiles_path):
    """Compute Morgan fingerprints for all 20 amino acids once."""
    with open(smiles_path) as f:
        smiles = json.load(f)

    lookup = {}
    for aa, smi in smiles.items():
        mol = Chem.MolFromSmiles(smi)
        fp = np.array(
            Chem.GetMorganFingerprintAsBitVect(mol, radius=FP_RADIUS, nBits=FP_NBITS)
        )
        lookup[aa] = fp

    zero_fp = np.zeros(FP_NBITS, dtype=int)
    lookup[' '] = zero_fp
    lookup['0'] = zero_fp
    return lookup


def compute_fingerprint(seq, aa_fps):
    """Compute the fingerprint feature vector for one sequence.

    Returns a 1D array of length MAX_LENGTH (mean-pooled across FP_NBITS).
    """
    fp = np.array([aa_fps[aa] for aa in seq])  # (seq_len, 2048)
    if len(seq) < MAX_LENGTH:
        padding = np.zeros((MAX_LENGTH - len(seq), FP_NBITS))
        fp = np.concatenate([fp, padding], axis=0)
    return np.mean(fp, axis=1)


# ---------------------------------------------------------------------------
# 3. Graph construction from ESM outputs
# ---------------------------------------------------------------------------

def build_graph_from_esm(rep, contact_map, seq_len):
    """Build a PyG Data object from ESM representation and contact map.

    Parameters
    ----------
    rep : torch.Tensor      (seq_len, 1280) per-token representation
    contact_map : np.ndarray (seq_len, seq_len) predicted contacts
    seq_len : int            original sequence length

    Returns
    -------
    DATA.Data  PyG graph object
    """
    # Pad representation to MAX_LENGTH.
    if seq_len < MAX_LENGTH:
        padding = torch.zeros(MAX_LENGTH - seq_len, ESM_DIM)
        feature = torch.cat([rep, padding], dim=0)
    else:
        feature = rep

    # Pad contact map to MAX_LENGTH x MAX_LENGTH, add self-loops.
    contact = np.zeros((MAX_LENGTH, MAX_LENGTH))
    contact[:seq_len, :seq_len] = contact_map
    contact = contact + np.eye(MAX_LENGTH)

    # Edge index.
    rows, cols = np.where(contact >= CONTACT_THRESHOLD)
    edge_index = torch.LongTensor(np.stack([rows, cols], axis=0))

    gcn_data = DATA.Data(
        x=feature.float(),
        edge_index=edge_index,
        y=torch.FloatTensor([0.0]),  # dummy label
    )
    gcn_data.__setitem__('target_size', torch.LongTensor([seq_len]))

    return gcn_data


# ---------------------------------------------------------------------------
# 4. Process one chunk: ESM -> graphs -> GNN predict -> results
# ---------------------------------------------------------------------------

class SimpleDataset(torch.utils.data.Dataset):
    def __init__(self, data_pro, data_fp):
        self.data_pro = data_pro
        self.data_fp = data_fp

    def __len__(self):
        return len(self.data_pro)

    def __getitem__(self, idx):
        return self.data_pro[idx], self.data_fp[idx]


def collate_fn(data_list):
    """Collate function matching the original collate()."""
    batchA = Batch.from_data_list([d[0] for d in data_list])
    batchB = default_collate([d[1] for d in data_list])
    return batchA, batchB


def process_chunk(
    chunk_ids, chunk_seqs, aa_fps,
    esm_model, alphabet, gnn_model, device,
    tokens_per_batch, gnn_batch_size,
):
    """Run the full pipeline on one chunk and return a DataFrame of results.

    Memory for graphs is freed before returning.
    """
    data_pro = []
    data_fp = []
    valid_seqs = []

    dataset = esm.FastaBatchedDataset(chunk_ids, chunk_seqs)
    batches = dataset.get_batch_indices(
        tokens_per_batch, extra_toks_per_seq=1
    )
    loader = torch.utils.data.DataLoader(
        dataset,
        collate_fn=alphabet.get_batch_converter(),
        batch_sampler=batches,
    )

    with torch.no_grad():
        for labels, strs, toks in loader:
            toks = toks.to(device, non_blocking=True)
            out = esm_model(
                toks, repr_layers=[ESM_LAYER], return_contacts=True
            )

            reps = out["representations"][ESM_LAYER].cpu()
            cts = out["contacts"].cpu()

            for i, label in enumerate(labels):
                seq_len = len(strs[i])
                seq = strs[i]

                rep = reps[i, 1 : seq_len + 1].clone()
                ct = cts[i, :seq_len, :seq_len].numpy()

                gcn_data = build_graph_from_esm(rep, ct, seq_len)
                fp = compute_fingerprint(seq, aa_fps)

                data_pro.append(gcn_data)
                data_fp.append(fp)
                valid_seqs.append(seq)

            del reps, cts, out

    if device.type == "cuda":
        torch.cuda.empty_cache()

    # -- GNN prediction on this chunk --
    ds = SimpleDataset(data_pro, data_fp)
    pred_loader = torch.utils.data.DataLoader(
        ds,
        batch_size=gnn_batch_size,
        shuffle=False,
        collate_fn=collate_fn,
    )

    gnn_model.eval()
    all_preds = []
    with torch.no_grad():
        for data_pro_batch, data_fp_batch in pred_loader:
            data_pro_batch = data_pro_batch.to(device)
            data_fp_batch = data_fp_batch.to(device)
            output = gnn_model(data_pro_batch, data_fp_batch)
            all_preds.append(output.cpu().numpy().flatten())

    mic_log10 = np.concatenate(all_preds)
    mic = 10.0 ** mic_log10

    # Free graph objects before returning.
    del data_pro, data_fp, ds, pred_loader
    if device.type == "cuda":
        torch.cuda.empty_cache()

    return pd.DataFrame({
        "sequence": valid_seqs,
        "MIC": mic,
        "MIC_unit": "uM",
    })


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = ArgumentParser(description="AMPredictor streaming inference")
    parser.add_argument("input", type=str, help="Path to input FASTA file")
    parser.add_argument("output", type=str, help="Path to output TSV file")
    parser.add_argument(
        "--tokens-per-batch", type=int, default=2048,
        help="Max tokens per ESM batch (lower = less GPU memory, default 2048)",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=5000,
        help="Sequences per processing chunk (lower = less CPU memory, default 5000)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=4096,
        help="Batch size for GNN prediction (default 4096)",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # --- 1. Read FASTA ---
    print(f"Reading FASTA from {args.input}...", file=sys.stderr)
    ids, seqs = read_fasta(args.input)
    print(
        f"  {len(ids)} sequences after filtering to <= {MAX_LENGTH} aa",
        file=sys.stderr,
    )

    if len(ids) == 0:
        pd.DataFrame(
            columns=["sequence", "MIC", "MIC_unit"]
        ).to_csv(args.output, sep="\t", index=False)
        print("No valid sequences. Empty output written.", file=sys.stderr)
        return

    # --- 2. Precompute fingerprint lookup ---
    smiles_path = os.path.join(script_dir, "data", "smiles_file.json")
    aa_fps = precompute_aa_fingerprints(smiles_path)

    # --- 3. Load both models once ---
    print("Loading ESM-1b model...", file=sys.stderr)
    esm_model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    esm_model.eval()
    esm_model = esm_model.to(device)

    print("Loading GNN model...", file=sys.stderr)
    model_path = os.path.join(script_dir, "models", "model_GNNPredictor_.model")
    gnn_model = GNNPredictor()
    gnn_model.load_state_dict(torch.load(model_path, map_location=device))
    gnn_model = gnn_model.to(device)

    # --- 4. Stream: chunk -> ESM -> graphs -> predict -> write -> free ---
    n_total = len(ids)
    header_written = False

    for chunk_start in range(0, n_total, args.chunk_size):
        chunk_end = min(chunk_start + args.chunk_size, n_total)
        chunk_ids = ids[chunk_start:chunk_end]
        chunk_seqs = seqs[chunk_start:chunk_end]

        chunk_df = process_chunk(
            chunk_ids, chunk_seqs, aa_fps,
            esm_model, alphabet, gnn_model, device,
            tokens_per_batch=args.tokens_per_batch,
            gnn_batch_size=args.batch_size,
        )

        chunk_df.to_csv(
            args.output,
            sep="\t",
            index=False,
            mode="a" if header_written else "w",
            header=not header_written,
        )
        header_written = True

        print(
            f"  Chunk {chunk_end}/{n_total} done "
            f"({chunk_end} sequences written)",
            file=sys.stderr,
        )

    del esm_model
    if device.type == "cuda":
        torch.cuda.empty_cache()

    print(
        f"Done: {n_total} predictions written to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()