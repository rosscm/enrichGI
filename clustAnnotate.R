# Annotate genes in pathways / complex
library(openxlsx)
source("readPathways.R")

######
# SET USER INPUTS
######

# Paths to gene-level qGI data + pathway annotation files
dataDir <- "/Users/catherineross/GIN"
inDir   <- sprintf("%s/bin/R/huGI/output_data/20200115", dataDir)

# qGI data
qGI_file <- list.files(pattern="PiScores", path=inDir, full.names=TRUE)
fdr_file <- list.files(pattern="FDR", path=inDir, full.names=TRUE)

# Pathway annotation file
path_file1 <- sprintf("%s/anno/Human_GOBP_AllPathways_no_GO_iea_October_01_2018_symbol.gmt", dataDir)
path_file2 <- sprintf("%s/anno/coreComplexes_Human_October_25_2018.tab", dataDir)

######
# START ANALYSIS
######

# Read in qGI data and filter
qGI <- read.delim(qGI_file, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
fdr <- read.delim(fdr_file, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

# Read in pathway file as a list
minGene=0
maxGene=500
path_file = path_file2

path <- readPathways(path_file, MIN_SIZE=minGene, MAX_SIZE=maxGene)

# Set up dataframes
neg <- data.frame(complex=NA, GI=NA, size=NA, n_size=NA, genes=NA)
pos <- data.frame(complex=NA, GI=NA, size=NA, n_size=NA, genes=NA)


##########
##### NOTE incomplete -- need to work on this #####
##########


# Fill in data
for (k in seq_along(path)) {
  path_test <- path[[k]]
  path_name <- names(path[k])

  dat_path_neg <- intersect(neg_dat, path_test)
  if (length(dat_path_neg) >= 2) {
    neg[k,]$complex <- path_name
    neg[k,]$GI <- "neg"
    neg[k,]$size <- length(path_test)
    neg[k,]$n_size <- length(dat_path_neg)
    neg[k,]$genes <- paste(dat_path_neg, collapse=", ")
  }

  dat_path_pos <- intersect(pos_dat, path_test)
  if (length(dat_path_pos) >= 2) {
    pos[k,]$complex <- path_name
    pos[k,]$GI <- "pos"
    pos[k,]$size <- length(path_test)
    pos[k,]$n_size <- length(dat_path_pos)
    pos[k,]$genes <- paste(dat_path_pos, collapse=", ")
  }
}

neg <- na.omit(neg)
pos <- na.omit(pos)
neg <- neg[order(neg$n_size, decreasing=TRUE),]
pos <- pos[order(pos$n_size, decreasing=TRUE),]
all <- list()
all[["negativeGI"]] <- neg
all[["positiveGI"]] <- pos
write.xlsx(all, "corumAnno_gi_PiScores_20200116_VPS52_241_standardCutoff.xlsx", row.names=FALSE)
