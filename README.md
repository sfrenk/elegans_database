# Get information about a list of genes/regions from Bionformatics databases

## biomart_lookup.R
R script to find features such as alternative IDs or Sequences for a list of genes. This script uses the Parasite (Wormbase) Biomart.

## samtools_lookup.py
Get sequence for regions in a bed file. Can also be loaded as a module to retrieve sequence of given coordinates. Requires local genome fasta file.

Examples:

Getting gene names and descriptions from a list of wormbase IDs
```
biomart_lookup.R -q wormbase_id -r gene_name,description -o results.txt input_ids.txt
```


Getting sequences from a bed file:
```
samtools_lookup.py -o seqs.fa coords.bed
```

Getting sequences from coordinates
```python
import samtools_lookup
my_seq = get_seq("I", 0, 101, name = None, zero = True, species = "elegans")

# Get sequence as a string
print(str(my_seq.seq))

```