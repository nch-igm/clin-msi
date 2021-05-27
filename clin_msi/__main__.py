import os
import logging
import argparse

from .train import train
from .predict import predict


def file_path_type(path: str) -> str:
    if not os.path.exists(path):
        raise ValueError(f"{path} does not exist.")
    return path


def comma_str_to_list(str_list):
    return str_list.split(",")


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='Workflow Type', dest="subparser_name")

    # train configuration
    train_parser = subparsers.add_parser("train")
    train_parser.add_argument('--input-bam-list', type=comma_str_to_list, required=True, help="input file with BAM paths and MSI status")

    # predict configuration
    predict_parser = subparsers.add_parser("predict")
    predict_parser.add_argument('--bam', type=str, required=True, help="path to BAM file")
    predict_parser.add_argument('--sample-name', type=str, required=True)
    predict_parser.add_argument('--model-dir', type=str, required=True, help="path to directory containing .pkl files")
    
    # shared configuration
    for subparser in [train_parser, predict_parser]:
        subparser.add_argument('--input-file', type=file_path_type, required=True, help="path to the tab separated input file with MSI regions")
        subparser.add_argument('--reference', type=file_path_type, required=True, help="path to reference(.fa or .fasta) file with index (.fai) in the same directory")
        subparser.add_argument('--output-dir', type=str, required=True, help="")
        subparser.add_argument('--allow-mismatch', action='store_true', help="allows a single base mismatch within the repeat region")
        subparser.add_argument('--normalization-scheme', type=str, required=True, choices=['z','std_u'], help="")

    args = parser.parse_args()

    if args.subparser_name == "train":
        train(
            args.input_file,
            args.input_bam_list,
            args.reference,
            args.allow_mismatch,
            args.normalization_scheme,
            args.output_dir)
    elif args.subparser_name == "predict":
        predict(
            args.input_file,
            args.bam,
            args.reference,
            args.allow_mismatch,
            args.sample_name,
            args.normalization_scheme,
            args.output_dir,
            args.model_dir)
    else:
        raise ValueError(f"Invalid workflow type {args.subparser_name}.")


if __name__ == "__main__":
    main()