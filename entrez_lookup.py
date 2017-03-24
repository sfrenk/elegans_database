#!/usr/bin/env python3

###############################################################################
# Get sequences for regions in a bed file
###############################################################################

import argparse
from Bio import Entrez, SeqIO
import os

###############################################################################
# Entrez setup
###############################################################################

# Entrez GI numbers for each chromosome:
chromosome_ids = {"I": 453231596, "II": 453231901, "III": 453232067, "IV": 453232348, "V": 453232767, "X": 453232919}

# Let Entrez know who you are
Entrez.email = "stephen.frenk@gmail.com"

###############################################################################

def get_seq(chrom, start, end, name = None, zero = False):
	
	if zero:
		# If entry is zero indexed (eg bed file), add 1 to the start coordinate
		start += 1

	if chrom in chromosome_ids.keys():

		handle = Entrez.efetch(db = "nucleotide", id = str(chromosome_ids.get(chrom)), rettype = "fasta", strand = 1, seq_start = start, seq_stop = end)
		record = SeqIO.read(handle, "fasta")
		handle.close()
		
		# Create fasta ID for entry
		if name != None:
			seq_id = ">" + name
		else:
			seq_id = ">" + chrom + ":" + str(start) + ".." + str(end)

		# Output results in fasta format
		return(seq_id + "\n" + str(record.seq))


def iterate_bed(input_file, output_file, zero = False):
	
	bed = open(input_file, "r")
	outfile = open(output_file, "w")
	
	# Iterate through bed file
	for line in bed:

		chrom = line.strip().split("\t")[0]
		# Remove "chr" prefix if it exists
		chrom = chrom.replace("chr", "")
		start = int(line.strip().split("\t")[1])
		end = int(line.strip().split("\t")[2])

		try:
			name = line.strip().split("\t")[3]
		except(IndexError):
			name = None

		seq = get_seq(chrom, start, end, name, zero = True)
		
		# Output results in multifasta format
		outfile.write(seq + "\n")

	bed.close()
	outfile.close()


def main():

	parser = argparse.ArgumentParser(description = 'Get sequences for list of regions in a bed file')

	parser.add_argument('file', help = 'file containing regions of interest in bed format')
	parser.add_argument('-o', '--output', help = 'output file name Default: entrez_results.txt', default = 'entrez_results.txt')

	args = parser.parse_args()

	iterate_bed(args.file, args.output, zero = True)


if __name__ == '__main__':
	main()
