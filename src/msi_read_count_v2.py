import pysam
import argparse
import re
from operator import itemgetter

def clin_msi_argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input-file', type=str, required=True)
    parser.add_argument('--bam', type=str, required=True)
    parser.add_argument('--reference', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True, help="")

    return parser

def repetitions(s):
   r = re.compile(r"(.+?)\1+")
   for match in r.finditer(s):
       yield (match.group(1), len(match.group(0))/len(match.group(1)))

def parse_input_file(input_file):
    location_list = []
    open_file = open(input_file, 'r')
    for line in open_file:
        chr, start, stop = line.split('\t')
        location_list.append([str(chr), int(start), int(stop)])

    return location_list

def clin_msi():
    parser = clin_msi_argparser()
    args = parser.parse_args()

    msi_location_list = parse_input_file(args.input_file)

    bam_file = pysam.AlignmentFile(args.bam, "rb")

    for chr, start, stop in msi_location_list:
        fasta_seq = pysam.faidx(args.reference, f"{chr}:{start-10}-{stop+10}").split('\n')[1]
        #fastq_seq = fastq_seq.split('\n')[1]
        rep_list = list(repetitions(fasta_seq))
        print(rep_list)
        largest_rep_unit = max(rep_list, key=itemgetter(1))[0]
        largest_rep_len = int(max(rep_list, key=itemgetter(1))[1])
        print(largest_rep_len, largest_rep_unit)
        search_string = fr"(\w{{5}}){largest_rep_unit}{{{largest_rep_len}}}(\w{{5}})"
        print(search_string)
        get_flanks = re.search(fr"(\w{{5}}){largest_rep_unit}{{{largest_rep_len}}}(\w{{5}})", fasta_seq)
        left_flank = get_flanks.group(1)
        right_flank = get_flanks.group(2)

        print(left_flank, fasta_seq, right_flank)


        for read in bam_file.fetch(chr, start-500, stop+500):
            if read.is_duplicate:
                continue
            if read.is_unmapped:
                continue
            if read.mapping_quality < 1:
                continue
            if read.reference_start > start:
                continue

            difference = start - read.reference_start

            #print(read.mapping_quality)
            #print(read)
            #print(read.query_sequence)
            read_start = read.query_alignment_start + difference
            #print(read.query_sequence[read_start:read.query_alignment_end])
            #print(list(repetitions(read.query_sequence[read_start:read.query_alignment_end])))
    #test = pysam.faidx("/Volumes/igm/apps/genomes/Homo_sapiens/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta", f"{chr}:{start}-{stop}")
    #print(list(repetitions(test.split('\n')[1])))

        #print(read.cigartuples)
        #if len(read.cigartuples) > 3:
            #print(read.cigartuples)
        #print(read.query_alignment_start)
        #print(read.query_alignment_end)
        #print(read.query_alignment_length)
        #print(read.reference_start)
        #print(read.cigartuples[0][0])
    #pysam.view(catch_stdout=False)
    #pysam.view("-q", "1", "-F", "1028", args.bam, "2:47641119-47641856")
    #iter = bam_file.fetch("2",47641419, 47641556)
    #for x in iter:
     #   print(str(x))


if __name__ == '__main__':
    clin_msi()