#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description = "Convert names: sequence, accession or cgc")

parser.add_argument("file", help = "input file (- for stdin)")
parser.add_argument("-o", "--output", help = "output file name (leave blank to print to stdout)", default = "-")
parser.add_argument("-q", "--query", help = "current name format")
parser.add_argument("-r", "--result", help = "format to convert to")

args = parser.parse_args()

if args.file == "-":
	args.file = sys.stdin

name_data = pd.read_csv("/home/sfrenk/Documents/Resources/Seq/WS251/xrefs.txt", sep = "\t", header = None,names = ["sequence", "accession", "cgc"])

infile = pd.read_csv(args.file, header = None, sep = "\t", squeeze = True)
#infile = pd.read_csv(args.file, header = None, sep = "\t")



###########
#for i in infile:
	#result =  name_data.loc[name_data["sequence"] == i, "cgc"]

#WARNING: the below command may crash the computer!
#outfile = infile.apply(lambda x: name_data.loc[name_data["sequence"] == x, "cgc"])
##########



outfile = name_data.loc[name_data[args.query].isin(infile), args.result]
#infile["converted_id"] = name_data.loc[name_data[args.query].isin(infile), args.result]

if args.output == "-":
	print(pd.DataFrame(outfile).to_string(header = False, index = False).replace( " ", ""))
	#print([row for row in outfile])

outfile.to_csv(args.output, sep = "\t", header = False, index = False)

#with open(args.output, "w") as f:
	#f.write([x + "\n" for x in outfile])