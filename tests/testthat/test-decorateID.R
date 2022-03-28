library(testthat)
library(Macarron)

prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
annotations_df = read.csv(file = prism_annotations, row.names = 1)
met_taxonomy = system.file("extdata", "demo_taxonomy.csv", package="Macarron")
taxonomy_df = read.csv(file = met_taxonomy)

# Checking output
decorate_output <- decorateID(annotations_df)
expect_equal(names(decorate_output), names(taxonomy_df))
expect_equal(dim(decorate_output), dim(taxonomy_df))
expect_gte(length(unique(decorate_output$Sub_Class)), length(unique(decorate_output$Class)))