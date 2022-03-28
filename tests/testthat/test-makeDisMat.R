library(testthat)
library(Macarron)

prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
abundances_df = read.csv(file = prism_abundances, row.names = 1)
prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
annotations_df = read.csv(file = prism_annotations, row.names = 1)
prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
metadata_df = read.csv(file = prism_metadata, row.names = 1)

# Checking min_prevalence
prism_mbx <- prepInput(abundances_df, annotations_df, metadata_df)
prism_w1 <- makeDisMat(prism_mbx) # default min prevalence
prism_w2 <- makeDisMat(prism_mbx, min_prevalence = 0.5) # lower min prevalence
expect_gt(dim(prism_w2)[2], dim(prism_w1)[2])