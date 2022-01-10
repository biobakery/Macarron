\name{MACARRoN}
\alias{MACARRoN}
\title{
MACARRoN (Metabolome Analysis and Combined Annotation Ranks for pRediction of
    Novel bioactives) is a workflow to systematically identify and prioritize 
    potentially bioactive (and often unannotated) small molecules in microbial 
    community metabolomic datasets. 
}
\description{
MACARRoN (Metabolome Analysis and Combined Annotation Ranks for pRediction of
    Novel bioactives) is a workflow to systematically identify and prioritize 
    potentially bioactive (and often unannotated) small molecules in microbial 
    community metabolomic datasets. It prioritizes metabolic features as 
    "potentially bioactive" using a combination of covariance with standard 
    i.e. annotated compounds, ecological properties such as prevalence and 
    relative abundance and association with environmental parameters or phenotypes.
}
\usage{
    MACARRoN(input_abundances, input_annotations, input_metadata, input_taxonomy)
}
\value{
    Results from applying the model.
}
\author{
    Amrisha Bhosle<abhosle@broadinstitute.org>,\cr
    Ludwig Geistlinger<Ludwig_Geistlinger@hms.harvard.edu>,\cr
    Maintainers: Sagun Maharjan<sagunmaharjann@gmail.com>,\cr
}
\examples{
    input_abundances <- system.file(
        'extdata','demo_abundances.csv', package="MACARRoN")
    input_annotations <-system.file(
        'extdata','demo_annotations.csv', package="MACARRoN")
    input_metadata <-system.file(
        'extdata','demo_metadata.csv', package="MACARRoN")
    input_taxonomy <-system.file(
        'extdata','demo_taxonomy.csv', package="MACARRoN")
    MACARRoN(input_abundances, input_annotations, input_metadata, input_taxonomy)
}