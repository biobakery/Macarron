library(testthat)
library(Macarron)

prism_abundances = system.file("extdata", "demo_abundances.csv", package="Macarron")
abundances_df = read.csv(file = prism_abundances, row.names = 1)
prism_annotations = system.file("extdata", "demo_annotations.csv", package="Macarron")
annotations_df = read.csv(file = prism_annotations, row.names = 1)
prism_metadata = system.file("extdata", "demo_metadata.csv", package="Macarron")
metadata_df = read.csv(file = prism_metadata, row.names = 1)
met_taxonomy = system.file("extdata", "demo_taxonomy.csv", package="Macarron")
taxonomy_df = read.csv(file = met_taxonomy)

prism_mbx <- prepInput(abundances_df, annotations_df, metadata_df)
prism_w <- makeDisMat(prism_mbx)

# Checking min_module_size
mms1 <- round(nrow(prism_w)^(1/3))
mms2 <- round(nrow(prism_w)^(1/3)) + 5
prism_modules <- findMacMod(prism_mbx, 
                            prism_w,
                            input_taxonomy = taxonomy_df,
                            min_module_size = mms1,
                            evaluateMOS = FALSE)
modules <- prism_modules[[1]]
n1 <- max(modules$module)
expect_gt(n1, 1)
prism_modules <- findMacMod(prism_mbx, prism_w,
                            input_taxonomy = taxonomy_df,
                            min_module_size = mms2,
                            evaluateMOS = FALSE)
modules <- prism_modules[[1]]
n2 <- max(modules$module)
expect_lt(n2, n1)