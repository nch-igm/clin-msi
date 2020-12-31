import re, os, sys, collections, glob, csv
from xlsxwriter.workbook import Workbook
from urllib.parse import urlparse
import argparse
import pysam

path_to_bam_file = sys.argv[1]
sample_name = sys.argv[2]
output_dir = sys.argv[3]
outputtable = open('{}/msi_count_table_{}.tsv'.format(output_dir, sample_name), 'w')

def msi_read_count_argparser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--input-file', action='append', nargs='+',choices={"male", "female"}, required=True)
    parser.add_argument('--bam', action='append', nargs='+', choices={"normal", "tumor"}, required=True)
    parser.add_argument('--reference', action='append', nargs='+', required=True)
    parser.add_argument('--output-dir', type=str, required=True, help="<REQUIRED> provide the s3 output path")

    return parser


def convert_tsv_to_xlsx(tsv_file):
    # Add some command-line logic to read the file names.
    xlsx_file = '{}/msi_count_table_{}.xlsx'.format(output_dir,sample_name)

    # Create an XlsxWriter workbook object and add a worksheet.
    workbook = Workbook(xlsx_file)
    worksheet = workbook.add_worksheet()

    # Create a TSV file reader.
    tsv_reader = csv.reader(open(tsv_file, 'r'), delimiter='\t')

    # Read the row data from the TSV file and write it to the XLSX file.
    for row, data in enumerate(tsv_reader):
        worksheet.write_row(row, 0, data)

    # Close the XLSX file.
    workbook.close()

def GoodfellowAnalysis(msiName, msiLocation, monomerLength, monomerType, monomerLocation, endCodon):
    bam_read_file = samtools(msiLocation, msiName)
    filtered_read_file = softclipping(bam_read_file, monomerLocation)
    monomerCount(monomerType, monomerLength, msiName, filtered_read_file, endCodon)


def samtools(msiLocation, msiName):
    bam_reads_dir = os.path.join(output_dir, 'bam_reads')
    read_dir = os.path.join(bam_reads_dir, msiName + 'READS')
    if not os.path.exists(read_dir): os.makedirs(read_dir)
    name = os.path.basename(path_to_bam_file).replace('.Recalibrated.bam', '_bam_reads')
    cmd = 'samtools view -q 1 -F 1028 %s %s > %s' %(path_to_bam_file, msiLocation, os.path.join(read_dir, name))
    os.system(cmd)
    return os.path.join(read_dir, name)

def softclipping(bam_read_file, monomerLocation):
    open_bam_read_file = open(bam_read_file, 'r')
    filteredoutput = open(bam_read_file.replace('_bam_reads', '_filtered_reads'), 'w')
    for line in open_bam_read_file:
        split_line = line.split('\t')
        difference = int(monomerLocation) - int(split_line[3])
        if int(split_line[3]) > int(monomerLocation): continue
        if re.search( r'^\d+S', split_line[5]):
            softclipstart = re.search( r'^(\d+)S', split_line[5])
            start = int(softclipstart.group(1)) + int(difference)
            if re.search( r'\d+S$', split_line[5]):
                softclipend = re.search( r'(\d+)S$', split_line[5])
                filteredoutput.write(split_line[9][start:-int(softclipend.group(1))] + '\n')
            else:
                filteredoutput.write(split_line[9][start:] + '\n')
        elif re.search( r'\d+S$', split_line[5]):
            softclipend = re.search( r'(\d+)S$', split_line[5])
            filteredoutput.write(split_line[9][difference:-int(softclipend.group(1))] + '\n')
        else:
            filteredoutput.write(split_line[9][difference:] + '\n')
    return bam_read_file.replace('_bam_reads', '_filtered_reads')

def monomerCount(monomerType, monomerLength, msiName, filtered_read_file, endCodon):
    finalTable = []
    numberSamples = 0
    total_counts = {}
    numberSamples += 1
    name = os.path.basename(filtered_read_file).replace('_filtered_reads', '')
    filtered_read_file = open(filtered_read_file, 'r')
    for line in filtered_read_file:
        if monomerType == 'A' and monomerLength >= 10:
            if line.endswith('AAAAAAA\n'): continue
            find_monomer = re.search( r'^[GCT]([A]+[TCG]?[A]{2,})%s' %(endCodon), line)
        elif monomerType == 'T':
            if line.endswith('TTTTTTT\n'): continue
            find_monomer = re.search( r'^[GCA]([T]+[ACG]?[T]{2,})%s' %(endCodon), line)
        elif monomerType == 'C':
            find_monomer = re.search( r'^[GAT]([C]+)%s' %(endCodon), line)
        elif monomerType == 'A' and monomerLength < 10:
            find_monomer = re.search( r'^[GCT]([A]+)%s' %(endCodon), line)
        elif monomerType == 'CT':
            find_monomer = re.search( r'^C((CT)\2*)%s' %(endCodon), line)
        else:
            print ('ERROR! Type not found!')
        if find_monomer:
            numNs = len(find_monomer.group(1))
            if monomerType == 'CT':
                numNs = int(len(find_monomer.group(1)) / 2)
            total_counts[numNs] = total_counts.get(numNs, 0) + 1

    sorted_total_counts = collections.OrderedDict(sorted(total_counts.items()))
    for key, value in sorted_total_counts.items():
        print (key, value)

    input_list = []
    input_list.append(name)
    for key in range(1, 31):
        if key in sorted_total_counts.keys():
            input_list.append(sorted_total_counts[key])
        else:
            input_list.append('0')

    finalTable.append(input_list)

    outputtable.write(msiName + '\n\n')
    for x in range(0, 31):
        for y in range(0, numberSamples):
            outputtable.write("%s\t" %(finalTable[y][x]))
        outputtable.write("\n")

if __name__ == '__main__':
    GoodfellowAnalysis('BAT26', '2:47641119-47641856', 27, 'A', 47641559, 'GGG')
    GoodfellowAnalysis('BAT25', '4:55597831-55598430', 25, 'T', 55598211, 'GAG')
    GoodfellowAnalysis('NR21', '14:23652229-23652649', 21, 'A', 23652346, 'GGC')
    GoodfellowAnalysis('NR24', '2:95849162-95849461', 23, 'T', 95849361, 'GTG')
    GoodfellowAnalysis('MONO27', '2:39536355-39536970', 27, 'T', 39536689, 'GAG')
    GoodfellowAnalysis('HSPH1_T17', '13:31722570-31722746', 17, 'A', 31722620, 'TCA')
    GoodfellowAnalysis('MSI01', '1:201754376-201754446', 17, 'T', 201754410, 'ATC')
    GoodfellowAnalysis('MSI03', '2:62063059-62063129', 17, 'A', 62063093, 'CTG')
    GoodfellowAnalysis('MSI04', '2:108479588-108479658', 18, 'T', 108479622, 'AA')
    GoodfellowAnalysis('MSI06', '5:172421726-172421796', 15, 'T', 172421760, 'AGA')
    GoodfellowAnalysis('MSI07', '6:142691916-142691986', 17, 'T', 142691950, 'AGC')
    GoodfellowAnalysis('MSI08', '7:1787485-1787555', 17, 'A', 1787519, 'GAC')
    GoodfellowAnalysis('MSI09', '7:74608706-74608776', 13, 'T', 74608740, 'ATG')
    GoodfellowAnalysis('MSI11', '11:106695477-106695550', 12, 'T', 106695514, 'GTA')
    GoodfellowAnalysis('MSI12', '15:45897737-45897807', 14, 'T', 45897771, 'CAG')
    GoodfellowAnalysis('MSI13', '16:18882625-18882695', 15, 'A', 18882659, 'CAG')
    GoodfellowAnalysis('MSI14', '17:19314883-19314953', 18, 'T', 19314917, 'CTC')
    GoodfellowAnalysis('MSH6', '2:48030589-48030689', 8, 'C', 48030639, 'TTC')
    GoodfellowAnalysis('ZFHX3', '16:72991718-72991777', 7, 'C', 72991757, 'TCC')
    GoodfellowAnalysis('CTCF', '16:67645312-67645371', 7, 'A', 67645338, 'CAA')
    GoodfellowAnalysis('ZFHX3_Dinucleotide', '16:72830894-72830953', 5, 'CT', 72830902, 'GGC')
    outputtable.close()
    convert_tsv_to_xlsx('{}/msi_count_table_{}.tsv'.format(output_dir, sample_name))




