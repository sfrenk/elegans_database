#!/usr/bin/env python

import os
from intermine.webservice import Service
import pandas as pd
import argparse
import sys

service = Service("http://im-dev.wormbase.org/tools/wormmine/service")

###############################################################################
# Python script for getting gene descriptions from Wormbase for a list of genes
# obtained from mmp pipeline
###############################################################################

# User varibles:

###############################################################################

parser = argparse.ArgumentParser(description='Convert between different nomenclature systems')

parser.add_argument('file', help='file contaning gene list or table')
parser.add_argument('-i', help='current format of gene list (wormbase, gene, transcript', required=True)
parser.add_argument('-o', help='format to convert to (wormbase, gene, transcript, description', required=True)
parser.add_argument('-c', type=int, default=0, help='index of column containing gene list (default: 0)')

args = parser.parse_args()

filename = args.file
gene_index = args.c

in_cmd = args.i
if in_cmd == 'wormbase':
	in_format = 'primaryIdentifier'
elif in_cmd == 'gene':
	in_format = 'symbol'
elif in_cmd == 'transcript':
	in_format = 'transcripts.symbol'
else:
	sys.exit('error: invalid input format')

out_cmd = args.o
if out_cmd == 'wormbase':
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

data = pd.read_csv(filename, sep='\t')

data[out_cmd+'_lookup'] = pd.Series(index=range(0, len(data)-1), dtype='str')

# Iterate through gene list and get brief description for each gene
for i in range(len(data)):
	gene = data.iloc[i, gene_index]
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
			data.ix[i, out_cmd+'_lookup'] = query_out
			print(query_out)    

data.to_csv(filename, sep='\t', index=False)
