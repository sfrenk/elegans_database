#!/usr/bin/env Rscript

library(biomaRt)
library(optparse)

option_list <- list(
    make_option(c("-q", "--query"),
                help = "format of query gene list (gene name, transcript name, etc)",
                type = "character",
                default = "gene_name"),
    make_option(c("-r", "--result"),
                help = "feature to look up (default = wormbase_id)",
                type = "character",
                default = "wormbase_id"),
    make_option(c("-c", "--column"),
                help = "column containing gene list (default = 1)",
                type = "integer",
                default = 1),
    make_option(c("-o", "--output"),
                help = "output file name (default = 'biomart_results.txt'",
                type = "character",
                default = "biomart_results.txt"),
    make_option(c("-e", "--header"),
                help = "use this flag if input file contains a header",
                type = "logical",
                action = "store_true",
                dest = "header",
                default = FALSE)
)

parser <- OptionParser(option_list = option_list, description = "Find features such as alternative IDs or Sequences for a list of genes")

arguments <- parse_args(parser, positional_arguments = 1)

opts <- arguments$options

in_term <- ifelse(
            opts$query == "gene", "gene_name", ifelse(
            opts$query == "wormbase_id", "ensembl_gene_id", ifelse(
            opts$query == "transcript", "ensembl_transcript_id", print("ERROR: invalid query value")
            )))

if (opts$result == "coords"){
    out_term <- c("chromosome_name", "start_position", "end_position", "strand")
} else {
    out_term <- ifelse(
                opts$result == "gene", "ensembl_gene_id", ifelse(
                opts$result == "wormbase_id", "ensembl_gene_id", ifelse(
                opts$result == "transcript", "wormbase_gseq", ifelse(
                opts$result == "type", "gene_biotype", ifelse(
                opts$result == "description", "description",
                print("ERROR: invalid result value")
                )))))
}

input_data <- read.table(arguments$args, sep = "\t", header = ifelse(opts$header, TRUE, FALSE))

mart <- useMart("parasite_mart", dataset = "wbps_eg_gene", host = "parasite.wormbase.org")

results <- getBM(attributes = c("ensembl_gene_id", out_term), 
                               filters = c("species_id_1010", in_term), 
                               values = c("celegans", input_data[opts$c]),
                               mart = mart)

input_data <- cbind(input_data, results)

write.table(input_data, file = opts$output, quote = FALSE, row.names = FALSE, col.names = ifelse(opts$header, TRUE, FALSE))
