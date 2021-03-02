#!/usr/local/bin/bash
# Script to run enrichGI on specified set of gene scores (eg crispRscore output)

# Define user parameters
## 1) input_type: one of chemgen (chemical), musgen (mouse), mgin (mouse GI), abgen (antibody), or pace
## 2) input_set: name of child folder where gene scores are stored
## 3) input_path: path to folder where gene scores are stored
## 4) anno_path: path to folder where annotation files are stored
input_type="gin"
input_set="20210107"
input_path="/Users/catherineross/data/qGI"
anno_path="/Users/catherineross/anno"

# I/O
dt=`date +%Y%m%d`
indir=${input_path}/${input_set}
outdir=output/${input_set}_${input_type}

# Specify appropriate gene set annotation file for CRISPR screen type
if [ $input_type = "musgen" ] || [ $input_type = "mgin" ]; then
  annotation_file=${anno_path}/MOUSE_GO_bp_no_GO_iea_symbol.gmt
  input_file=${indir}/*_masterTable.txt
else
  #annotation_file=${anno_path}/Human_GO_bp_no_GO_iea_October_01_2018_symbol.gmt
  annotation_file=${anno_path}/coreComplexes_Corum_Human_October_25_2018.gmt
  input_file=${indir}/qGI_${input_set}.txt
fi

# Generate dated output directory
if [ -d "$outdir" ]; then
  echo "Output directory '$outdir' exists"
else
  mkdir $outdir
fi

# Run enrichGI on given inputs
Rscript enrichGI.R -i ${input_file} \
                   -a ${annotation_file} \
                   --MIN_GENE 0 \
                   --MAX_GENE 500 \
                   -o ${outdir} \
                   -q 149,307,335
