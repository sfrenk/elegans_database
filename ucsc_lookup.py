#!/usr/bin/env python3

from urllib.request import urlopen
import xml.etree.ElementTree as ET
import os
import re
import pandas as pd

# Example url for getting c elegans genome sequence data from ucsc: http://genome.ucsc.edu/cgi-bin/das/ce10/dna?segment=II:100000,200000

###############################################################################
# Python script for getting gene descriptions from Wormbase for a list of genes
# obtained from mmp pipeline
###############################################################################

# User varibles:

# Working directory:
wd = '/media/sfrenk/Acer/Users/Stephen/Documents/Sequencing/mmp/gwas'

# Name of file containing top hits

hitsfile = 'rdna_top_0.05__pval.txt'

# Number of top hits to query (based on poisson p-values)
top = 30

# output file name
output = 'sequences.txt'

###############################################################################
print('Looking up sequences for region list:')


data = pd.read_csv(hitsfile, sep='\t', nrows=top)

sequences = []

for i in range(len(data)):
	chromosome = data.loc[i, 'chromosome']
	region = data.loc[i, 'region']
	print(region)
	start = str(int(region) - 500)
	end = str(int(region) + 499)
	sequence = urlopen('http://genome.ucsc.edu/cgi-bin/das/ce10/dna?segment='+chromosome+':'+start+','+end).read()
	# This gives an xml file - the 'DNA' element needs to be extracted from this
	sequence = ET.fromstring(sequence)
	sequence = sequence.findall('./SEQUENCE/DNA')
	# Line breaks will be included as '\n' in the sequence so we need to get rid of those
	sequence = re.sub('\\\\n', '', str(sequence))
	sequences.append(chromosome+':'+start+'-'+end+':	'+str(sequence))

f = open(output, 'w')

for line in sequences:
	f.write(line+'\n')

f.close()

print('results written to '+output)
