
"""
        How to use
------------------------------
Python version = 3.7

python parseRaw.py example_raw_data.xlsx

example_raw_data.xlsx = path to the raw output generated by Patrick's Scripts

Should take < 1 minutes to run

"""

import argparse
import pandas as pd
import numpy as np
import os
import re

def parse_raw_data(path, sample_name):

	# Fixed microsatellite markers
	markers = ['BAT26','BAT25','NR21','NR24','MONO27','HSPH1_T17','MSI01','MSI03','MSI04','MSI06','MSI07','MSI08','MSI09','MSI11','MSI12','MSI13','MSI14','MSH6','ZFHX3','CTCF','ZFHX3_Dinucleotide']
	markers_1 = ['BAT25','NR21','NR24','MONO27','HSPH1_T17','MSI01','MSI03','MSI04','MSI06','MSI07','MSI08','MSI09','MSI11','MSI12','MSI13','MSI14','MSH6','ZFHX3','CTCF','ZFHX3_Dinucleotide', 'end']
	marker_zip = list(zip(markers, markers_1))
	marker_list = []
	
	# Generate the input for the design matrix
	output_intermediate = '/'.join(path.split('/')[:-1])+'/intermediate.csv'
	try:
		raw_df = pd.read_excel(path)
	except pd.errors.ParserError:
		print('Input is not in excel format trying csv') 
		try:
			raw_df = pd.read_csv(path)
			if 'Unnamed: 0' in raw_df.columns:
				raw_df = raw_df.drop('Unnamed: 0',axis=1)
		except pd.errors.ParserError:
			print('Input is not in csv. Parsing Failed')
		finally:
			return

	rows_to_drop = np.where(pd.notnull(raw_df).sum(axis=1) == 0)[0]
	raw_df = raw_df.drop(rows_to_drop)
	end_capped_row = raw_df.index.max() + 1 # add a token that acts like an EOF 
	raw_df.loc[end_capped_row] = ['end_capped'] + list(np.repeat(np.nan,raw_df.shape[1]-1))
	raw_df.columns = [raw_df.columns[0]] + ['' for i in range(len(raw_df.columns[1:]))]
	raw_df.to_csv(output_intermediate, index=False)

	# data must be read as a text file to parse
	content = []
	with open(output_intermediate, 'r') as fin:
		for line in fin:
			line = line.strip()
			content.append(line)

	# continue to parse by extracting each marker
	content = "\r".join(content)
	for m in marker_zip:
		marker_list.append(re.search("{}.*?(?={})".format(m[0],m[1]),content).group())

	# if using python2 we might need to use this code
	# with open(dilution_table, 'r') as fin:
	# 	for line in fin: # only 1 line
	# 		line = line.strip()
	# 		for m in marker_zip:
	# 			marker_list.append(re.search("{}.*?(?={})".format(m[0],m[1]),line).group())
	# 		# split into markers
			
	# build the dataframe
	list_of_df = []
	for m in marker_list:
		feat = m.split("\r")[:-1] # extract the relevant fields
		mrk = feat[0].split(",")[0] # extract the maker id
		samples = feat[1].split(",") # extract the samples
		value_list = [] # store the reads at each nucleotide
		for f in feat[2:]:
			value_list.append(f.split(","))
		# convert list into a dataframe. (reads are in order of sampels)	
		value_list = np.asarray(value_list).T
		z_score_val = []
		# zscore normalize the attributes by row
		for val in value_list:
			val = np.asarray(val).astype(int)
			mu = np.mean(val)
			std = np.std(val)
			if val.sum() != 0:
				val = val - mu / std
				z_score_val.append(val)
			else:
				z_score_val.append(val)
		z_score_val = np.asarray(z_score_val)
		columns = [str(i) + "_" +mrk for i in range(0,value_list.shape[1])]
		list_of_df.append(pd.DataFrame(data=z_score_val, index=samples, columns=columns))

	# Join all the matrices
	marker_design_matrix = list_of_df[0]
	for df in list_of_df[1:]:
		marker_design_matrix = marker_design_matrix.join(df)

	# save the matrix
	final_output = '/'.join(path.split('/')[:-1])+'/zscore_normalized_{}.csv'.format(sample_name)
	marker_design_matrix.to_csv(final_output)

def main():
	parser = argparse.ArgumentParser(prog='Parse raw MSI data',description='Convert an excel file into a design matrix for machine learning with the markers zscore normalized by row')
	parser.add_argument('path', metavar='/path/to/zscore/normalized/msi/data/sample.xlsx', type=str, nargs=1, help='path to the test data must be in excel or csv format')
	parser.add_argument('sample_name', metavar='XXXX', type=str, nargs=1, help='name of sample', default='test_sample')
	args = parser.parse_args()
	path = args.path[0]
	sample_name = args.sample_name[0]
	parse_raw_data(path, sample_name)

if __name__ == '__main__':
    main()
