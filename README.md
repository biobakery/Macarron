

# MACARRoN User Manual #  
  
MACARRoN (Metabolome Analysis and Combined Annotation Ranks for pRediction of Novel bioactives) is a workflow to systematically identify and prioritize potentially bioactive (and often unannotated)
small molecules in microbial community metabolomic datasets. It prioritizes metabolic features as "potentially bioactive" using a combination of covariance with standard i.e. annotated compounds, ecological properties such as prevalence and relative abundance and association with environmental parameters or phenotypes.

If you use the MACARRoN software, please cite our manuscript: 

If you have questions, please direct it to:
[MACARRoN Forum](https://forum.biobakery.org/c/microbial-community-profiling/macarron)

--------------------------------------------
   
## Contents ##
* [Description](#description)
* [Requirements](#requirements)
* [Installation](#installation)

## Description

MACARRoN (Metabolome Analysis and Combined Annotation Ranks for pRediction of Novel bioactives) is a workflow to systematically identify and prioritize potentially bioactive (and often unannotated) small molecules in microbial community metabolomic datasets.

## Requirements

MACARRoN is an R package that can be run as a set of R functions or via a command line wrapper.
     
## Installation

If only running from the command line, you do not need to install the MACARRoN package but you will need to install the MACARRoN dependencies.

### From the command line ###

1. Download the source: [MACARRoN.master.zip] (https://github.com/biobakery/MACARRoN/archive/master.zip)
2. Decompress the download:
    * ``$ tar xzvf MACARRoN-master.zip``
3. Install the Bioconductor dependencies: SummarizedExperiment, BiocParallel, DelayedArray, WGCNA, Maaslin2 and ComplexHeatmap.
4. Install the CRAN dependencies:
    * ``$ R -q -e "install.packages(c('optparse','logging','data.table','ff','ffbase','dynamicTreeCut','plyr','psych','ggplot2','circlize','grDevices'), repos='http://cran.r-project.org')"``
5. Install the MACARRoN package (only required if running as an R function).

### From R ###

To install the latest version of MACARRoN:

```{r, eval=FALSE}
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MACARRoN")
```
To install the latest development version of MACARRoN:

```{r, eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("biobakery/MACARRoN")
```

