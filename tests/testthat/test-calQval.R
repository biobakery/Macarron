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
prism_modules <- findMacMod(prism_mbx, 
                            prism_w,
                            input_taxonomy = taxonomy_df,
                            evaluateMOS = FALSE)
prism_ava <- calAVA(prism_mbx,
                    prism_modules)

# Checking output
prism_qval <- calQval(prism_mbx,
                      prism_modules,
                      plot_heatmap = FALSE)
expect_equal(length(unique(prism_qval$metadata)), 1)
expect_equal(length(unique(prism_qval$value)), length(unique(prism_mbx[["diagnosis"]])) - 1)
