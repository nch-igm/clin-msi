import pysam
import argparse
import os
import regex as re
from operator import itemgetter
from collections import defaultdict
import pandas as pd
from .normalize_counts import parse_raw_data
from .apply_msi_model import apply_model

def clin_msi_argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input-file', type=str, required=True, help="path to the tab separated input file with MSI regions")
    parser.add_argument('--bam', type=str, required=True, help="path to BAM file")
    parser.add_argument('--reference', type=str, required=True, help="path to reference(.fa or .fasta) file with index (.fai) in the same directory")
    parser.add_argument('--sample-name', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True, help="")
    parser.add_argument('--model-dir', type=str, required=True, help="path to directory containing .pkl files")
    parser.add_argument('--allow-mismatch', action='store_true', help="allows a single base mismatch within the repeat region")
    parser.add_argument('--normalization-scheme', type=str, required=True, choices=['z', 'std_u'], help="")

    return parser

def repeat_finder(s):
    #Taken from https://stackoverflow.com/questions/9079797/detect-repetitions-in-string
    r = re.compile(r"(.+?)\1+")
    for match in r.finditer(s, overlapped=True):
        yield (match.group(1), len(match.group(0))/len(match.group(1)))

def parse_input_file(input_file):
    location_list = []
    open_file = open(input_file, 'r')
    for line in open_file:
        chr, start, stop = line.split('\t')
        location_list.append([str(chr), int(start), int(stop)])

    return location_list

def predict():
    parser = clin_msi_argparser()
    args = parser.parse_args()

    msi_location_list = parse_input_file(args.input_file)

    #Load BAM file
    bam_file = pysam.AlignmentFile(args.bam, "rb")
    df = pd.DataFrame()

    for chr, start, stop in msi_location_list:
        #get expected repeat length from provided reference file
        length_dict = defaultdict(int)
        fasta_seq = pysam.faidx(args.reference, f"{chr}:{start-10}-{stop+10}").split('\n')[1]
        rep_list = list(repeat_finder(fasta_seq))
        largest_rep_unit = max(rep_list, key=itemgetter(1))[0]
        largest_rep_len = int(max(rep_list, key=itemgetter(1))[1])
        get_flanks = re.search(fr"(\w{{5}})(?:{largest_rep_unit}){{{largest_rep_len}}}(\w{{5}})", fasta_seq)
        left_flank = get_flanks.group(1)
        right_flank = get_flanks.group(2)

        #Skip reads that are unmapped, duplicate, or qual=0
        for read in bam_file.fetch(chr, start, stop):
            if read.is_duplicate:
                continue
            if read.is_unmapped:
                continue
            if read.mapping_quality < 1:
                continue
            if read.reference_start > start:
                continue

            read_wo_softclip = read.query_sequence[read.query_alignment_start:read.query_alignment_end]

            #allow one mismatched base in repeat region if specified
            if args.allow_mismatch:
                get_rep = re.search(fr"{left_flank}({largest_rep_unit}+[ACTG]?{largest_rep_unit}+){right_flank}", read_wo_softclip)
            else:
                get_rep = re.search(fr"{left_flank}({largest_rep_unit}+){right_flank}", read_wo_softclip)

            # if repeat is not found, go to next read
            if get_rep is None:
                continue

            length_dict[len(get_rep.group(1))] += 1

        length_list = ['N-10', 'N-9', 'N-8', 'N-7', 'N-6', 'N-5', 'N-4', 'N-3', 'N-2', 'N-1', 'N',
                       'N+1', 'N+2', 'N+3', 'N+4', 'N+5', 'N+6', 'N+7', 'N+8', 'N+9', 'N+10']
        repeat_count_list = []
        for x in range(largest_rep_len - 10, largest_rep_len + 11):
            if not length_dict[x]:
                length_dict[x] = 0
            repeat_count_list.append(length_dict[x])

        df['Repeat_Length'] = length_list
        df[f'{chr}:{start}-{stop}'] = repeat_count_list
    df.to_csv('test_raw_count.csv', index=False)
    normalized_df = parse_raw_data(df, args.sample_name, args.normalization_scheme)

    normalized_df.to_csv(os.path.join(args.output_dir, args.sample_name + '_normalized.csv'), index=False)

    #apply model to normalized msi counts
    final_results_file = os.path.join(args.output_dir, args.sample_name + '_MSIscore.txt')
    apply_model(os.path.join(args.output_dir, args.sample_name + '_normalized.csv'), args.model_dir, final_results_file)

if __name__ == '__main__':
    predict()