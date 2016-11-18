#!/usr/bin/env python3

###############################################################################
# Get sequences for regions in a bed file
###############################################################################

import argparse
from Bio import Entrez, SeqIO
import os

###############################################################################
# ARGUMENTS
###############################################################################

parser = argparse.ArgumentParser(description = 'Get sequences for list of regions in a bed file')

parser.add_argument('file', help = 'file containing regions of interest')
parser.add_argument('-o', '--output', help = 'output file name Default: entrez_results.txt', default = 'entrez_results.txt')
#parser.add_argument('-c', '--coords', help = 'query database with chromosomal coordinates. Coordinates should be in bed format (tab separated list of chromosome, start and end for each line)', action = 'store_true', default = False)

args = parser.parse_args()

in_file = open(args.file, "r")

if os.path.exists(args.output):
	os.remove(args.output)
out_file = open(args.output, "a")

# Entrez GI numbers for each chromosome:
chromosome_ids = {"I": 453231596, "II": 453231901, "III": 453232067, "IV": 453232348, "V": 453232767, "X": 453232919}

# Let Entrez know who you are
Entrez.email = "stephen.frenk@gmail.com"

###############################################################################
###############################################################################

# Iterate through bed file
for line in in_file:

	chrom = line.strip().split("\t")[0]
	# Remove "chr" prefix if it exists
	chrom = chrom.replace("chr", "")
	start = line.strip().split("\t")[1]
	end = line.strip().split("\t")[2]

	# Ignore header if it exists
	if chrom in chromosome_ids.keys():

		print(line.strip())

		handle = Entrez.efetch(db = "nucleotide", id = str(chromosome_ids.get(chrom)), rettype = "fasta", strand = 1, seq_start = start, seq_stop = end)
		record = SeqIO.read(handle, "fasta")
		handle.close()
		
		# Output results in multifasta format
		out_file.write(">" + chrom + ":" + start + ".." + end + "\n" + str(record.seq) + "\n")

in_file.close()
out_file.close()
