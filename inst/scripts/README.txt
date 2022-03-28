++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Running Macarron in command line
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Authors: Amrisha Bhosle, Ludwig Geistlinger, Sagun Maharjan

Contents:
  1. Usage
  2. Demo example

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Usage
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Maccaron requires four input files. Make sure to provide the full path to the Macarron executable.

$ Rscript MacarronCMD.R --help

Output: 
Usage: MacarronCMD.R [Inputs]
 <feature_abundances.csv>
 <feature_annotations.csv>
 <sample_metadata.csv>
 <chemical_taxonomy.csv> 

Options:
	-h, --help
		Show this help message and exit

	-o OUTPUT, --output=OUTPUT
		Output folder name [Default: Macarron_output]

	-m METADATA_VARIABLE, --metadata_variable=METADATA_VARIABLE
		Metadata column with bioactivity-relevant epidemiological or environmental conditions (factors) 
                [Default: Second column of metadata CSV file]

	-p MIN_PREVALENCE, --min_prevalence=MIN_PREVALENCE
		Minimum percent of samples of each condition in <metadata variable>
                in which metabolic feature is recorded [Default: 0.7]

	-e EXECUTION_MODE, --execution_mode=EXECUTION_MODE
		BiocParallel Class for execution [Default: serial] [Choices:serial, multi]

	-k STANDARD_IDENTIFIER, --standard_identifier=STANDARD_IDENTIFIER
		Annotation such as HMDB/METLIN/PubChem ID for an annotated metabolite feature 
                [Default: Second column of annotation CSV file]

	-a ANCHOR_ANNOTATION, --anchor_annotation=ANCHOR_ANNOTATION
		Metabolite name for an annotated metabolite feature 
                [Default: Third column of annotation CSV file]

	-c MIN_MODULE_SIZE, --min_module_size=MIN_MODULE_SIZE
		Minimum module size for module generation 
                [Default: cube root of number of prevalent features]

	-f FIXED_EFFECTS, --fixed_effects=FIXED_EFFECTS
		Maaslin2 parameter; Fixed effects for differential abundance linear model, 
                comma-delimited for multiple effects [Default: all]

	-r RANDOM_EFFECTS, --random_effects=RANDOM_EFFECTS
		Maaslin2 parameter; Random effects for differential abundance linear model, 
                comma-delimited for multiple effects [Default: none]

	-d REFERENCE, --reference=REFERENCE
		Maaslin2 parameter; The factor to use as a reference for a variable with more than two levels 
                provided as a string of 'variable,reference' semi-colom delimited for multiple variables 
                [Default: Alphabetically first]

	-n CORES, --cores=CORES
		Maaslin2 parameter; The number of R processes to run in parallel [Default: 1]

	-l PLOT_HEATMAP, --plot_heatmap=PLOT_HEATMAP
		Maaslin2 parameter; Generate a heatmap for significant associations

	-s PLOT_SCATTER, --plot_scatter=PLOT_SCATTER
		Maaslin2 parameter; Generate scatter plots for significant associations

	-t HEATMAP_FIRST_N, --heatmap_first_n=HEATMAP_FIRST_N
		Maaslin2 parameter; Generate heatmap for top n significant associations [Default: 50]

	-b SHOW_BEST, --show_best=SHOW_BEST
		Write 1000 or fewer highly prioritized metabolic features 
                 into a separate file [Default: TRUE]

	-w PRIORITY_THRESHOLD, --priority_threshold=PRIORITY_THRESHOLD
		Cut-off of priority score for showing highly
                 prioritized features [Default: 0.9]

	-x PER_MODULE, --per_module=PER_MODULE
		Show first n highly prioritized features in a module [Default: 10]

	-y PER_PHENOTYPE, --per_phenotype=PER_PHENOTYPE
		Show highly prioritized n features per phenotype/condition [Default: 1000]

	-z ONLY_CHARACTERIZABLE, --only_characterizable=ONLY_CHARACTERIZABLE
		Show highly prioritized features in modules which contain at least
                 one annotated metabolite [Default: TRUE]

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Demo example
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Example input files can be found in the inst/extdata folder of Macarron source. To run the demo, make sure to give the full path to the executable and the demo input files:

$ Rscript MacarronCMD.R demo_abundances.csv \
demo_annotations.csv \
demo_metadata.csv \
demo_taxonomy.csv

