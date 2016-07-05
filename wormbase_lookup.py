#!/usr/bin/env python3

import os
from intermine.webservice import Service
import pandas as pd
import argparse
import sys

service = Service("http://intermine.wormbase.org/tools/wormmine/service")

###############################################################################
# Python script for getting gene descriptions/IDs from Wormbase for a list of genes
###############################################################################

# User varibles:

###############################################################################

parser = argparse.ArgumentParser(description='Look up useful fetures for a gene list in wormbase and add these as a new column')

parser.add_argument('file', help='file contaning gene list or table')
parser.add_argument('-q', help='query: current format of gene list (wormbase_id, gene, transcript). Default: transcript', default='transcript')
parser.add_argument('-r', help='result: feature to look up (wormbase_id, gene, transcript, description). Default: description', default='description')
parser.add_argument('-c', type=int, default=0, help='index of column containing gene list. Default: 0')
parser.add_argument('--header', help='input file contains a header', action='store_true')
parser.add_argument('-o', help='output file name Default: wormbase_results.txt', default='wormbase_results.txt')

args = parser.parse_args()

filename = args.file
gene_index = args.c

in_cmd = args.q
if in_cmd == 'wormbase_id':
	in_format = 'primaryIdentifier'
elif in_cmd == 'gene':
	in_format = 'symbol'
elif in_cmd == 'transcript':
	in_format = 'transcripts.symbol'
else:
	sys.exit('error: invalid input format')

out_cmd = args.r
if out_cmd == 'wormbase_id':
	out_format = 'primaryIdentifier'
elif out_cmd == 'gene':
	out_format = 'symbol'
elif out_cmd == 'transcript':
	out_format = 'transcripts.symbol'
elif out_cmd == 'description':
	out_format = 'briefDescription'
elif out_cmd == 'length':
	out_format = 'sequence.length'
else:
	sys.exit('error: invalid output format')

if args.header:
	head = 0
else:
	head = None

outfile = args.o

data = pd.read_csv(filename, sep='\t', header=head)

data[out_cmd] = pd.Series(index=range(0, len(data)-1), dtype='str')

# Rename query column to 'query_in'
new_cols = data.columns.values
new_cols[0] = 'query_in'
data.columns = new_cols

# Iterate through gene list and get brief description for each gene
for i in range(len(data)):
	gene = data.iloc[i, gene_index]
	gene = gene.strip()
	query = service.new_query("Gene")
	# The view specifies the output columns
	query.add_view(out_format)

	try:
		query.add_constraint(in_format, "=", gene, code = "A")
	# If the primary identifier doesn't work, try a different identifier
	except intermine.model.ModelError:
		try:
			query.add_constraint('secondaryIdentifier', "=", gene, code = "A")
		except intermine.model.ModelError:
			print('could not convert: ', gene)
		else: query_success = True
	else:
		query_success = True
		
	if query_success:	
		for row in query.rows():
			query_out = (row[out_format])
			data.ix[i, out_cmd] = query_out
			print(gene, " ", query_out)

data = data[['query_in', out_cmd]]

data.to_csv(outfile, sep='\t', index=False)
