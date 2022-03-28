library(testthat)
library(Macarron)

prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
abundances_df = read.csv(file = prism_abundances, row.names = 1)
prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
annotations_df = read.csv(file = prism_annotations, row.names = 1)
prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
metadata_df = read.csv(file = prism_metadata, row.names = 1)

# Samples lacking either feature abundance or metadata information
metadata_df_sub <- head(metadata_df, 75) # picking only 75 samples
prism_mbx <- prepInput(abundances_df, annotations_df, metadata_df_sub)
expect_lt(dim(SummarizedExperiment::assay(prism_mbx))[2],dim(abundances_df)[2])

abundances_df_sub <- abundances_df[,1:75] # picking only 75 samples
prism_mbx <- prepInput(abundances_df_sub, annotations_df, metadata_df)
expect_lt(dim(SummarizedExperiment::assay(prism_mbx))[2],dim(abundances_df)[2])