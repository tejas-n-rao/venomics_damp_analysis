import argparse

import benchmark_utils as utils


def main():
    parser = argparse.ArgumentParser(description="Run inference for MBC-Attention")
    parser.add_argument("-i", "--input_file", dest="input_file", type=str,
                        help="Path to the input FASTA file.", required=True)
    parser.add_argument("-o", "--output_file", dest="output_file", type=str,
                        help="Path to the inference output TSV file.", required=True)
    parser.add_argument("-m", "--model_file", dest="model_file", type=str,
                        help="Path to the MBC-Attention model file.", required=False,
                        default="model/whole_train.mdl")

    args = parser.parse_args()

    processed_inputs, sequence_df = utils.convert_inputs(args.input_file)
    pred_df = utils.run_inference(
        processed_inputs, sequence_df, args.model_file
    )
    utils.convert_outputs(pred_df, args.output_file)

if __name__ == "__main__":
    main()
