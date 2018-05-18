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

# Fasta files for query genomes

# LOCAL
genome_fasta_files = {"elegans" : "/home/sfrenk/Documents/Resources/Seq/WS251/genome/genome.fa", "briggsae" : "/home/sfrenk/Documents/Resources/Seq/briggsae/wormbase/WS263/genome.fa", "remani" : "/home/sfrenk/Documents/Resources/Seq/remani/genome.fa"}

# CLUSTER
#genome_fasta_files = {"elegans" : "/nas/longleaf/home/sfrenk/proj/seq/WS251/genome/genome.fa", "briggsae" : "/proj/ahmedlab/steve/seq/briggsae/WS263/caenorhabditis_briggsae.PRJNA10731.WBPS10.genomic.fa", "remani" : "/proj/ahmedlab/steve/seq/briggsae/WS263/caenorhabditis_remanei.PRJNA53967.WBPS10.genomic.fa.gz"}

###############################################################################
def get_seq(chrom, start, end, name = None, zero = False, species = "elegans"):
	
	# Get coords
	if zero:
		# If entry is zero indexed (eg bed file), add 1 to the start coordinate
		start = int(start) + 1

	coords = str(chrom) + ":" + str(start) + "-" + str(end)

	# Get genome fasta for species

	# Get sequence
	print("Getting sequence from " + str(species) + " at " + str(chrom) + ":" + str(start) + "-" + str(end) + "\n")
	genome_fasta = genome_fasta_files[species]
	seq = subprocess.check_output(["samtools", "faidx", genome_fasta, coords]).decode("utf-8") 
		
	if name == None:
		name = str(re.search(">([^\n]+)\n", seq).group(1))

	# Remove fasta header
	seq = re.sub(">[^\n]*\n", "", seq)
	seq = re.sub("\n", "", seq).upper()

	# Output result as seqrecord object
	seq_record = SeqRecord(Seq(seq, generic_dna), id = name, description = "")

	return(seq_record)


def iterate_bed(input_file, output_file, zero = True, species = "elegans"):

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

		output_seq = get_seq(chrom, start, end, name, zero = True, species = species)

		# Output results in multifasta format
	
		outfile.write(output_seq.format("fasta") + "\n")


	bed.close()
	outfile.close()


def main():

	parser = argparse.ArgumentParser(description = 'Get sequences for list of regions in a bed file')

	parser.add_argument('file', help = 'file containing regions of interest in bed format')
	parser.add_argument('-o', '--output', help = 'output file name Default: samtools_results.txt', default = 'samtools_results.txt')
	parser.add_argument('-s', '--species', help = 'elegans (default), briggsae or remani', default = 'elegans')

	args = parser.parse_args()

	iterate_bed(args.file, args.output, species = args.species)


if __name__ == '__main__':
	main()
