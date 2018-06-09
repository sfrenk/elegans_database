#!/usr/bin/env Rscript
# version 3.1

suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(optparse))

connect_to_biomart <- function(biomart_site = "ensembl", biomart_dataset = "celegans_gene_ensembl"){
    # Paraiste isn't working for me, using Esembl instead
    #mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "parasite.wormbase.org")
    # There's sometimes a problem with connecting to biomart, so need to keep trying if there is an error
    try_counter <- 1
    while(!exists("mart_object") & try_counter < 10){
        print(paste0("Connection attempt ", as.character(try_counter), "..."))
        tryCatch(mart_object <- useMart(biomart = biomart_site, dataset = biomart_dataset), error = function(x) {Sys.sleep(10)})
        try_counter <- try_counter + 1
    }
    # Exit if it still hasnt't worked after 10 attempts
    if (!exists("mart_object")){
        print("ERROR: Could not connect to Biomart")
        return(0)
    } else {
        print("Connected!")
        return(mart_object)
    }
}

get_biomart_results <- function(filter_term, atts, input_list, biomart_object){
    
    # The filter term is added to the output attributes to allow merging of results with query list.
    # This allows fors multiple attribute results from a single query and makes sure the results match up with the input
    atts <- c(filter_term, atts)
    
    # Get results from biomart
    results <- getBM(attributes = atts, 
                     filters = filter_term, 
                     values = input_list,
                     mart = biomart_object)
    
    # Merge filter and attributes. This ensures that the results remain in the same order as the query
    input_genes <- data.frame(input_list)
    colnames(input_genes) <- as.character(atts[1])
    input_genes$index_col <- 1:nrow(input_genes)
    
    results_table <- merge(input_genes, results, by = atts[1], all.x = TRUE)
    results_table <- results_table[order(results_table$index_col, decreasing = FALSE),]
    results_table$index_col <- NULL
    results_table <- results_table[!(duplicated(results_table[atts[1]])),]
    
    # The "description" attribute comes with source info, which isn't very useful, so we can get rid of this.
    if ("description" %in% colnames(results_table)){
        suppressPackageStartupMessages(library(stringr))
        results_table$description <- sapply(results_table$description, function(x) str_replace(x, '[ ]*\\[Source:.*', ""))
    }
    return(results_table)
}

output_data <- function(results_table, output){
    
    # Output data
    if (output != "STDOUT"){
        write.table(results_table, file = output, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        print(results_table, right = FALSE, row.names = FALSE)
    }
}

main <- function(){
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
                    type = "character",
                    default = "STDOUT"),
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
                opts$query == "gene", "external_gene_name", ifelse(
                opts$query == "wormbase_id", "wormbase_gene", ifelse(
                opts$query == "seq_id", "wormbase_gene_seq_name", ifelse(
                opts$query == "transcript", "ensembl_transcript_id", print("ERROR: invalid query value")
                ))))

    opts$result <- unlist(strsplit(opts$result, ","))

    
    if ("gene" %in% opts$result){
        out_term <- c(out_term, "external_gene_name")
    }
    if ("wormbase_id" %in% opts$results){
        out_term <- c(out_term, "wormbase_gene") 
    }
    if ("seq_id" %in% opts$results){
        out_term <- c(out_term, "wormbase_gene_seq_name") 
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
    if ("go" %in% opts$results){
        out_term <- c(out_term, "GO term name")
    }

    write("LOOKING UP:", stderr())
    write(out_term, stderr())
    write("FROM:", stderr())
    write(in_term, stderr())

    # Read data from standin if "-" option specified, otherwise read from file
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
    
    # Try to connect to biomart, exit if this fails
    mart <- connect_to_biomart()
    
    if (!exists("mart")){
        print("Exiting due to biomart connection failure")
        quit(save = "no", status = 1, runLast = FALSE)
    }
    
    # Get biomart results
    results <- get_biomart_results(filters = in_term, atts = out_term, input_list = input_data[opts$c])
    
    # Output data
    output_data(results_table = results, output = opts$output)
    
    print("biomart_lookup: Done!")
        
}

