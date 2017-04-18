#!/usr/bin/env python3

###############################################################################
# Get sequences for regions in a bed file using samtools faidx
# This is a local alternative to entrez_lookup that's much faster
###############################################################################

import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import os
import subprocess
import re

###############################################################################
# VARIABLES

# Fasta file of query genome
genome_fasta = "/home/sfrenk/Documents/Resources/Seq/WS251/genome/genome.fa"

###############################################################################

def get_seq(chrom, start, end, name = None, zero = False):
	
	if zero:
		# If entry is zero indexed (eg bed file), add 1 to the start coordinate
		start += 1

	coords = str(chrom) + ":" + str(start) + "-" + str(end)
	seq = subprocess.check_output(["samtools", "faidx", genome_fasta, coords]).decode("utf-8") 
		
	if name == None:
		name = str(re.search(">([^\n]+)\n", seq).group(1))

	# Remove fasta header

	seq = re.sub(">[^\n]*\n", "", seq)
	seq = re.sub("\n", "", seq).upper()

	# Output result as seqrecord object
	seq_record = SeqRecord(Seq(seq, generic_dna), id = name, description = "")

	return(seq_record)


def iterate_bed(input_file, output_file, zero = True):

	""" Get sequences for regions in a bed file
	"""
	
	bed = open(input_file, "r")

	# Initialize output file
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

		output_seq = get_seq(chrom, start, end, name, zero = True)

		# Output results in multifasta format
	
		outfile.write(output_seq.format("fasta") + "\n")


	bed.close()
	outfile.close()


def main():

	parser = argparse.ArgumentParser(description = 'Get sequences for list of regions in a bed file')

	parser.add_argument('file', help = 'file containing regions of interest in bed format')
	parser.add_argument('-o', '--output', help = 'output file name Default: samtools_results.txt', default = 'samtools_results.txt')

	args = parser.parse_args()

	iterate_bed(args.file, args.output)


if __name__ == '__main__':
	main()
