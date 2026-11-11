import argparse

import benchmark_utils as utils


def main():
    parser = argparse.ArgumentParser(description="Run inference for Deep-AMP")
    parser.add_argument("-i", "--input_file", dest="input_file", type=str,
                        help="Path to the input FASTA file.", required=True)
    parser.add_argument("-o", "--output_file", dest="output_file", type=str,
                        help="Path to the inference output TSV file.", required=True)
    parser.add_argument("-m", "--model_file", dest="model_file", type=str,
                        help="Path to the Deep-AMP model file.", required=True,
                        choices=[
                            "saved_models/paper/regressor/CNN_gr_neg",
                            "saved_models/paper/regressor/CNN_gr_pos",
                            "saved_models/paper/regressor/LSTM_gr_neg",
                            "saved_models/paper/regressor/LSTM_gr_pos"
                        ])

    args = parser.parse_args()
    sequence_mapping, sequences = utils.convert_inputs(args.input_file)
    mic = utils.run_inference(sequences, args.model_file)
    utils.convert_outputs(sequence_mapping, mic, args.output_file)

if __name__ == "__main__":
    main()