library(testthat)
source("../biomart_lookup.R")

test_that("biomart_actually_connects", {
    
    # Make sure biomart actually connects
    mart <- connect_to_biomart()
    expect_that(typeof(mart), equals("S4"))

    # Check that results are as expected
    results <- get_biomart_results(filter_term = c("external_gene_name"), atts = c("wormbase_gene"), input_list = c("trt-1"), biomart_object = mart)
    expect_that(as.character(results["wormbase_gene"]), equals("WBGene00006618"))

})

