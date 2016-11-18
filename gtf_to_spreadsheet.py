import re
import os

gtf_file = "/home/sfrenk/Documents/Resources/Seq/WS251/genes.gtf"
gene_file = "/home/sfrenk/Documents/Database/genes.txt"
transcript_file = "/home/sfrenk/Documents/Database/transcripts.txt"

if os.path.exists(gene_file):
	os.remove(gene_file)
if os.path.exists(transcript_file):
	os.remove(transcript_file)

gtf = open(gtf_file, "r")
genes = open(gene_file, "a")
transcripts = open(transcript_file, "a")

for line in gtf:
	print(line)
	atr = line.strip().split("\t")
	chromosome = atr[0]
	start = atr[3]
	end = atr[4]
	strand = atr[6]
	gene_id = re.search("gene_id ([^ ;]*)(;|$)", atr[8]).group(1)
	gene_name = re.search("gene_name ([^ ;]*)(;|$)", atr[8]).group(1)
	wormbase_accession = re.search("wormbase_accession ([^ ;]*)(;|$)", atr[8]).group(1)
	transcript_id = re.search("transcript_id ([^ ;]*)(;|$)", atr[8]).group(1)
	
	genes.write("\t".join([gene_id, gene_name, wormbase_accession, chromosome, start, end, strand, "\n"]))
	transcripts.write("\t".join([transcript_id, gene_id, "\n"]))


gtf.close()
genes.close()
transcripts.close()
