# Macarron User Manual #  
  
Macarron (Metabolome Analysis and Combined Annotation Ranks for pRediction of Novel bioactives) is a workflow to systematically identify and prioritize potentially bioactive (and often unannotated)
small molecules in microbial community metabolomic datasets. It prioritizes metabolic features as "potentially bioactive" using a combination of covariance with standard i.e. annotated compounds, ecological properties such as prevalence and relative abundance and association with environmental parameters or phenotypes.

If you use the Macarron software, please cite our manuscript: 

If you have questions, please direct it to:
[Macarron Forum](https://forum.biobakery.org/c/microbial-community-profiling/Macarron)

--------------------------------------------
   
## Contents ##
* [Description](#description)
* [Requirements](#requirements)
* [Installation](#installation)

## Description

Macarron (Metabolome Analysis and Combined Annotation Ranks for pRediction of Novel bioactives) is a workflow to systematically identify and prioritize potentially bioactive (and often unannotated) small molecules in microbial community metabolomic datasets.

## Requirements

Macarron is an R package that can be run as a set of R functions or via a command line wrapper.
     
## Installation

If only running from the command line, you do not need to install the Macarron package but you will need to install the Macarron dependencies.

### From the command line ###
1. Install the Bioconductor dependencies: SummarizedExperiment, BiocParallel, DelayedArray, Maaslin2 and ComplexHeatmap.
2. Install the CRAN dependencies:
    * ``$ R -q -e "install.packages(c('optparse','logging','data.table','WGCNA','ff','ffbase','dynamicTreeCut','plyr','psych','ggplot2','circlize','grDevices','xml2','RCurl','RJSONIO'), repos='http://cran.r-project.org')"``

### From R ###

To install the latest version of Macarron:

```{r, eval=FALSE}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Macarron")
```
To install the latest development version of Macarron:

```{r, eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("biobakery/Macarron")
```

