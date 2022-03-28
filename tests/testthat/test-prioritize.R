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
prism_qval <- calQval(prism_mbx,
                      prism_modules)
prism_es <- calES(prism_mbx,
                  prism_qval)

# Checking output
prism_prioritized <- prioritize(prism_mbx,
                                prism_modules,
                                prism_ava,
                                prism_qval,
                                prism_es)
expect_true(is.data.frame(prism_prioritized[[1]]))
expect_true(is.data.frame(prism_prioritized[[2]]))
expect_equal(ncol(prism_prioritized[[1]]), ncol(prism_prioritized[[2]]))
expect_equal(ncol(prism_prioritized[[1]]), ncol(annotations_df) + 10)
             