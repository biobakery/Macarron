# Macarron User Manual #  
  
Macarron (Metabolome Analysis and Combined Annotation Ranks to pRioritize Novel bioactives) is a workflow to systematically identify and prioritize potentially bioactive (and often unannotated)
small molecules in microbial community metabolomic datasets. Macarron prioritizes metabolic features as potentially bioactive in a phenotype/condition of interest using a combination of (a) covariance with annotated metabolites, (b) ecological properties such as abundance with respect to covarying annotated compounds, and (c) differential abundance in the phenotype/condition of interest.

If you have questions, please direct it to:
[Macarron Forum](https://forum.biobakery.org/c/microbial-community-profiling/Macarron)

## Installation ##
Macarron requires R version 4.2.0 or higher. Install Bioconductor and then install Macarron:

```
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Macarron")
```

## How to Run ##
Macarron can be run from the command line or as an R function. Both methods require the same
arguments, have the same options, and use the same default settings. The package includes the
wrapper `Macarron()` as well as functions which perform different steps in the Macarron
framework.

### Input CSV files ###
Macarron requires 4 comma-separated, appropriately formatted input files. The files and their 
formatting constraints are described below.

1. Metabolic features abundances
    * Must contain features in rows and samples in columns.
    * First column must identify features.
2. Metabolic features annotations
    * Must contain features in rows and annotations in columns.
    * First column must identify features.
    * Second column must contain either HMDB ID or PubChem Compound Identifier (CID).
    * Third column must contain the name of the metabolite.
    * Fourth column must contain a continuous chemical property such as m/z or RT or shift/ppm.
    * Other annotations such as RT, m/z or other identifiers can be listed column 4 onward.
3. Sample metadata
    * Must contain samples in rows and metadata in columns.
    * First column must identify samples.
    * Second column must contain categorical metadata relevant to prioritization such as phenotypes, exposures or environments.
4. Chemical taxonomy
    * First column must contain the HMDB ID or PubChem CID. IDs must be consistent between annotation and taxonomy files.
    * Second and third columns must contain chemical subclass and class of the respective metabolite.
    
If you do not have the chemical taxonomy file, you can generate this file using the annotation dataframe and Macarron utility `decorate_ID` (see Advanced Topics).

### Output Files ###
By default, all files will be stored in a folder named Macarron_output inside the current working directory. The main prioritization results are stored in ``prioritized_metabolites_all.csv``. Another file, ``prioritized_metabolites_characterizable.csv`` is a subset of ``prioritized_metabolites_all.csv`` and only contains metabolic features which covary with at least one annotated metabolite.
The columns in these output files are:

- Feature_index: Lists the identifier of the metabolic feature found in column 1 of abundance and annotation files.
- HMDB_ID (or PubChem ID): Public database identifier from column 2 of annotation file (column 1 of annotation dataframe).
- Metabolite name: From column 2 of annotation dataframe.
- mz: The continuous numerical chemical property from column 3 of the annotation dataframe.
- Priority_score: 1 indicates most prioritized. It is the percentile from the meta-rank of AVA, q-value and effect size.
- Status: Direction of perturbation (differential abundance) in the phenotype (or environment) of interest compared to reference phenotype.
- Module: ID of the covariance module a metabolic feature is a member of. Module = 0 indicates a singleton i.e., a metabolic feature that is not assigned to any module.
- Anchor (of a module): Metabolic feature that has the highest abundance in any phenotype. 
- Related_classes: Chemical taxonomy of the annotated features that covary with a metabolic feature.
- Covaries_with_standard: 1 (yes) and 0 (no). Column specifies if the metabolic feature covaries with at least one annotated (standard) metabolite.
- AVA: Abundance versus anchor which is a ratio of the highest abundance (in any phenotype) of a metabolic feature and highest abundance of the covarying anchor. Naturally, the AVA of an anchor metabolite is 1.
- qvalue: Estimated from multivariate linear model using `Maaslin2`.
- effect_size
- Remaining columns from the annotation dataframe are appended.

### Run a demo in R ###

#### Using CSV files as inputs ####
Example (demo) input files can be found under ``inst/extdata`` folder of the `Macarron` source. These files were generated from the [PRISM](https://pubmed.ncbi.nlm.nih.gov/30531976/) study of stool metabolomes of individuals with inflammatory bowel disease (IBD) and healthy "Control" individuals. Control and IBD are the two phenotypes in this example. Macarron will be applied to prioritize metabolic features with respect to their bioactivity in IBD. Therefore, in this example, the phenotype of interest is "IBD" and the reference phenotype is "Control". The four input files are ``demo_abundances.csv``, ``demo_annotations.csv``, ``demo_metadata.csv``, and ``demo_taxonomy.csv``. 

```
library(Macarron)
prism_abundances <- system.file(
    'extdata','demo_abundances.csv', package="Macarron")
prism_annotations <-system.file(
    'extdata','demo_annotations.csv', package="Macarron")
prism_metadata <-system.file(
    'extdata','demo_metadata.csv', package="Macarron")
mets_taxonomy <-system.file(
    'extdata','demo_taxonomy.csv', package="Macarron")
prism_prioritized <- Macarron::Macarron(input_abundances = prism_abundances,
                                        input_annotations = prism_annotations,
                                        input_metadata = prism_metadata,
                                        input_taxonomy = mets_taxonomy)
```

#### Using dataframes as inputs ####
```
abundances_df = read.csv(file = prism_abundances, row.names = 1) # setting features as rownames
annotations_df = read.csv(file = prism_annotations, row.names = 1) # setting features as rownames
metadata_df = read.csv(file = prism_metadata, row.names = 1) # setting samples as rownames 
taxonomy_df = read.csv(file = mets_taxonomy)

# Running Macarron
prism_prioritized <- Macarron::Macarron(input_abundances = abundances_df,
                                        input_annotations = annotations_df,
                                        input_metadata = metadata_df,
                                        input_taxonomy = taxonomy_df)
```

