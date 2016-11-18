#!/usr/bin/env python3

import subprocess
import os
import re
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = "convert wormbase gff3 annotations file to a gtf file that can be used in RNA-seq analysis")

parser.add_argument("file", help = "gff3 (input) filename")
parser.add_argument("-o", "--output", help = "gtf (output) filename")
parser.add_argument("-x", "--xrefs", help = "xrefs file")

args = parser.parse_args()

# The xrefs file contains name conversions for various different naming conventions
print("Loading xrefs data")
xrefs = pd.read_csv(args.xrefs, sep = "\t", comment = "/", header = None)
xrefs.columns = ["sequence", "accession", "gene", "transcript", "protein", "INSDC_parent", "INSDC_locus", "INSDC_protein", "uniprot"]

# Filter gff3 file for WormBase exon annotations
#print("filtering input file")
#with open("temp.wormbase", "w") as f:
#	subprocess.call(["grep", "WormBase", args.file], stdout = f)

#with open("temp.exon", "w") as f:
#	subprocess.call(["grep", "exon", "temp.wormbase"], stdout = f)

# Process gff3 file
print("Preparing to process input file")
if os.path.exists("genes.gtf"):
	os.remove("genes.gtf")

gff_file = open("temp.exon", "r")
gtf_file = open(args.output, "a")
excluded = open("gtf_excluded.txt", "w")


for line in gff_file:
	skip_line = False
	transcript_id = re.search("Parent=[^:]*:([^;]*)", line.strip().split()[8]).group(1)
	try:
		gene_id = str(xrefs.loc[xrefs.transcript == transcript_id, "sequence"].values[0])
	except(IndexError):
		try:
			transcript_id = re.search("Parent=[^:]*:([^.;]*[.]?[^.]*)", line.strip().split()[8]).group(1)
			gene_id = str(xrefs.loc[xrefs.transcript == transcript_id, "sequence"].values[0])
		except(IndexError):
			excluded.write(line)
			skip_line = True

	if not skip_line:
		gene_name = str(xrefs.loc[xrefs.transcript == transcript_id, "gene"].values[0])
		if gene_name == ".":
			gene_name = gene_id
		wormbase_accession = str(xrefs.loc[xrefs.transcript == transcript_id, "accession"].values[0])
		attribute = "gene_id " + gene_id + "; gene_name " + gene_name + "; transcript_id " + transcript_id + "; wormbase_accession " + wormbase_accession
		gtf_entry = line.strip().split()[0:8]
		gtf_entry.append(str(attribute))
		print(gtf_entry)
		gtf_file.write("\t".join(gtf_entry) + "\n")
	

gff_file.close()
gtf_file.close()
excluded.close()

os.remove("temp.wormbase")
