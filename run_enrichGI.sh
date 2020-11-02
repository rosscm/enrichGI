#!/usr/local/bin/bash
# Script to run enrichGI on specified set of gene scores (eg crispRscore output)

# Define user parameters
## 1) crispr_type: one of chemgen (chemical), musgen (mouse), mgin (mouse GI), abgen (antibody), or pace
## 2) crispr_path: path to folder where gene scores are stored
## 2) crispr_path: name of child folder where gene scores are stored
crispr_type="gin"
crispr_path="/Users/catherineross/projects/qGI/output_data"
#crispr_path="/Users/catherineross/projects/GIN/scripts/crispRscore/output"
crispr_set="20201028"

# I/O
dt=`date +%Y%m%d`
indir=${crispr_path}/${crispr_set}
outdir=output/${crispr_set}

# Specify appropriate gene set annotation file for CRISPR screen type
if [ $crispr_type = "musgen" ] || [ $crispr_type = "mgin" ]; then
  #annotation_file=input/Mouse_Human_Reactome_April_01_2019_symbol.gmt
  annotation_file=input/Mouse_Human_KEGG_December_24_2015_symbol.gmt
  input_file=${indir}/*_masterTable.txt
else
  #annotation_file=input/Human_Reactome_October_01_2018_symbol.gmt
  annotation_file=input/Human_KEGG_December_24_2015_symbol.gmt
  input_file=${indir}/qGI_${crispr_set}.txt
fi

# Generate dated output directory
if [ -d "$outdir" ]; then
  echo "Output directory '$outdir' exists"
else
  mkdir $outdir
fi

# Run enrichGI on score file
Rscript enrichGI.R -i ${input_file} \
                   -a ${annotation_file} \
                   --MIN_GENE 0 \
                   --MAX_GENE 10000 \
                   -o ${outdir} \
                   -q 287
