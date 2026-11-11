from Bio import SeqIO
from collections import OrderedDict
import numpy as np
import pandas as pd
from tensorflow import keras
from typing import List, Literal, Tuple

from code.utils import chopping, onehot_encoding, padding


def convert_inputs(input_file: str, alphabet: str = "BCDSQKIPTFNGHLRWAVEYM-"
) -> Tuple[OrderedDict, np.array]:
    """Convert inputs for the Deep-AMP regressor.

    Args:
        input_file: path to the input fasta file.
        alphabet: amino acid alphabet for one hot encoding. It should include
            begin and end tokens for padding.

    Returns:
        A tuple (sequence_mapping, sequences) where sequence_mapping is a
        subset of proteins retrieved from the input file that satisfied
        model's input requirements; and sequences is an array consisting
        of corresponding preprocessed protein sequences.
    """
    print("Converting inputs")
    sequence_to_id = OrderedDict()
    for sequence in SeqIO.parse(input_file, "fasta"):
        sequence_to_id[str(sequence.seq)] = sequence.id

    sequences = chopping(list(sequence_to_id.keys()))
    sequence_mapping = OrderedDict(
        {seq: sequence_to_id[seq] for seq in sequences}
    )

    sequences = padding(sequences)
    sequences = onehot_encoding(sequences, alphabet)
    sequences = np.asarray(sequences)

    print(f"Finished converting inputs. Discarded "
          f"{len(sequence_to_id.keys()) - len(sequence_mapping.keys())} sequences")
    return sequence_mapping, sequences


def run_inference(sequences: np.array, model_file: Literal[
    "saved_models/paper/regressor/CNN_gr_neg",
    "saved_models/paper/regressor/CNN_gr_pos",
    "saved_models/paper/regressor/LSTM_gr_neg",
    "saved_models/paper/regressor/LSTM_gr_pos"
]) -> List:
    """Run inference for the Deep-AMP regressor.

    Args:
        sequences: preprocessed protein sequences.
        model_file: path to one of the Deep-AMP regressors.

    Returns:
        A list of predicted log10 MIC values.
    """
    print("Running inference")
    model = keras.models.load_model(model_file, compile=False)
    preds = model.predict(sequences, batch_size=64).flatten().tolist()
    print("Finished inference")
    return preds

def convert_outputs(sequence_mapping: OrderedDict, mic: List, output_file: str):
    """Convert outputs of the Deep-AMP regressor.

    Args:
        sequence_mapping: sequence-id mapping of the proteins, which were not
            discarded during the convertion of outputs, and for which
            predictions were computed.
        mic: log10_MIC values predicted during inference.
        output_file: path to the output tsv file.
    """
    print("Converting outputs")
    df = dict()
    for (seq, seq_id), y in zip(sequence_mapping.items(), mic):
        df[seq_id] = {
            "Sequence": seq, "Prediction": "AMP" if y < 4 else "non-AMP",
            "Log10_MIC": y, "MIC": 10 ** y, "MIC_unit": "uM"
        }

    df = pd.DataFrame.from_dict(df, orient="index").reset_index(names="Sequence_id")
    df.to_csv(output_file, sep="\t", index=False)
    print("Finished converting outputs")
