## enrichGI

Pipeline for running gene-set enrichment analysis (GSEA) on genetic interaction data.

## Prerequisites

R version > 3.3

CRAN packages: data.table, dplyr, reshape2, ggplot2, scales, forcats, openxlsx, argparse

Bioconductor packages: fgsea

## Installation

After ensuring that the above R packages are installed, download this repository locally to any location.

## Usage

Navigate to this repo's folder in your command line, e.g. with `cd enrichGI`. Run this pipeline with
the command `Rscript enrichGI.R`, followed by the required arguments `-i` for the input qGI score file and `-a` for the pathway annotation file. See below for a description of how
the inputs should be formatted (see files provided in `input` for examples).

### Input

**Genetic interaction scores** (`-i`, `--input_file`): Text file with rows as GI scores per gene and columns as screens. *REQUIRED*.

**Pathway annotations** (`-a`, `--annotation_file`): GMT file with one row per gene set annotation. The file provided in `input` was downloaded from the [Bader Lab](http://baderlab.org/GeneSets) and includes all human Reactome pathway annotations.
*REQUIRED*: default "input/Human_Reactome_October_01_2018_symbol.gmt".

**Query** (`-q`, `--query`): Name or ID number of screen or screens to run GSEA
	(e.g., TAZ (will test all TAZ queries in data), TAZ_208, or 208). Input is *not* case sensitive.
*OPTIONAL*: default "all".

**Output folder** (`-o`, `--output_folder`): Path to an output folder for all results.
*OPTIONAL*: default "output".

**Min pathway size** (`--MIN_GENE`): Integer specifying minimum number of genes to be considered in a pathway.
*OPTIONAL*: default 10.

**Max pathway size** (`--MAX_GENE`): Integer specifying maximum number of genes to be considered in a pathway.
*OPTIONAL*: default 500.

**Enrichment FDR** (`--SIG_FDR`): Integer specifying FDR threshold to define significant pathway enrichment.
*OPTIONAL*: default 0.05.

*NO NEED TO WORRY ABOUT THE FOLLOWING PARAMETERS*

**Permutations** (`--SET_PERM`): Integer specifying number of permutations to run to determine enrichment significance.
	Minimal possible nominal p-value ~= 1/`SET_PERM`.
*OPTIONAL*: default "10000".

**Seed** (`--SET_SEED`): Integer specifying seed value to maintain result consistency. *OPTIONAL*: default "42".


### Example command
`Rscript enrichGI.R -i input/20200323/gi_PiScores_20200323_v2018_20200323.txt -q TAZ`
