#!/usr/bin/env Rscript
# version 2.1

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
                help = "output results to specified file name (by default, results are printed to terminal)",
                type = "character"),
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

# The structure below allows for multple "return" results
# Having the query filter in the first column of the results table allows for merging of dataframes later on

out_term <- ifelse(
    opts$query == "gene", "external_gene_id", ifelse(
    opts$query == "transcript", "ensembl_transcript_id", in_term))

opts$result <- unlist(strsplit(opts$result, ","))

if ("gene" %in% opts$result){
    out_term <- c(out_term, "external_gene_id")
}
if ("wormbase_id" %in% opts$results){
    out_term <- c(out_term, "ensembl_gene_id") 
}
if ("transcript" %in% opts$result){
    out_term <- c(out_term, "ensembl_transcript_id")
}
if ("type" %in% opts$result){
    out_term <- c(out_term, "gene_biotype")
}
if ("description" %in% opts$result){
    out_term <- c(out_term, "description")
}
if ("coords" %in% opts$result){
    out_term <- c(out_term, "chromosome_name", "start_position", "end_position", "strand")
}

print("LOOKING UP:")
print(out_term)
print("FROM:")
print(in_term)

if (arguments$args == "-") {
    input_data <- character(length = 0)
    f <- file("stdin")
    open(f)
    while(length(line <- readLines(f, n = 1)) > 0) {
        input_data <- c(input_data, line)
    }
    input_data <- data.frame("V1" = input_data)

} else {
    input_data <- read.table(arguments$args, sep = "\t", header = ifelse(opts$header, TRUE, FALSE))
}

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "parasite.wormbase.org")

results <- getBM(attributes = c(out_term), 
                               filters = c(in_term), 
                               values = c(input_data[opts$c]),
                               mart = mart)

# Merge filter and attributes. This ensures that the results remain in the same order as the query
input_genes <- data.frame(input_data[opts$c])
colnames(input_genes) <- as.character(out_term[1])

results_table <- merge(input_genes, results, by = out_term[1], sort = FALSE, all.x = TRUE)

if (length(opts$output) > 0){

    write.table(results_table, file = opts$output, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    print(results_table)
}