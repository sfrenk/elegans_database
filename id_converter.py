#!/usr/bin/env python3

# be sure to use Python 3!
import pandas as pd
import argparse

###############################################################################
# Python script for getting gene descriptions from Wormbase for a list of genes
# obtained from mmp pipeline
###############################################################################

# User varibles:
refs_file = "/media/sfrenk/Acer/Users/Stephen/Documents/Resources/Seq/WS251/refs.txt"

###############################################################################

parser = argparse.ArgumentParser(description='Look up useful fetures for a gene list in wormbase and add these as a new column')

parser.add_argument('file', help='file contaning gene list or table')
parser.add_argument('-q', help='query: current format of gene list (wormbase_id, gene, transcript). Default: transcript', default='transcript')
parser.add_argument('-r', help='result: feature to look up (wormbase_id, gene, transcript, description). Default: description', default='description')
parser.add_argument('-c', type=int, default=0, help='index of column containing gene list. Default: 0')
parser.add_argument('--header', help='input file contains a header', action='store_true')
parser.add_argument('-o', help='output file name Default: wormbase_results.txt', default='converter_results.txt')
parser.add_argument('--inplace', help='edit file in place rather than creating a new file', action='store_true' )

args = parser.parse_args()

filename = args.file
gene_index = args.c

in_cmd = args.q
if in_cmd == 'wormbase_id':
	in_format = 1
elif in_cmd == 'gene':
	in_format = 2
elif in_cmd == 'transcript':
	in_format = 3
else:
	sys.exit('error: invalid input format')

out_cmd = args.r
if out_cmd == 'wormbase_id':
	out_format = 1
elif out_cmd == 'gene':
	out_format = 2
elif out_cmd == 'transcript':
	out_format = 3
else:
	sys.exit('error: invalid output format')

if args.header:
	head = 0
else:
	head = None

outfile = args.o

if args.inplace:
	outfile = filename

data = pd.read_csv(filename, sep='\t', header=head)
refs = pd.read_csv("/media/sfrenk/Acer/Users/Stephen/Documents/Resources/Seq/WS251/refs.txt", comment='/', header=None, sep='\t')

for i in range(len(data)):
	gene = data.iloc[i, gene_index]
	hit = refs.loc[refs[in_format] == gene, [out_format]]
	try:
		hit = str(hit.iloc[0,0])
	except IndexError:
		hit = str(hit)
	# If there is no gene name for the feature, use default name
	if len(hit) <= 1:
		hit = refs.loc[refs[in_format] == gene, [0]]
		try:
			hit = str(hit.iloc[0,0])
		except IndexError:
			hit = str(hit)
	# If no hit can be found, keep the original entry
	# Note: often if no hit can be found, the hit returned will be "Empty DataFrame". This result can easily be detected due to its character length
	if len(hit) <= 1 or len(hit) > 37:
		hit = str(gene)
	data.iloc[i, gene_index] = hit
	print(gene, "---->", hit)

data.to_csv(outfile, sep='\t', index=False)
