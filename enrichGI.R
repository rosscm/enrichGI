#!/usr/bin/env Rscript

######
# LOAD PACKAGES
######

# Load packages
packages <- c("fgsea", "data.table", "dplyr", "reshape2", "ggplot2", "scales",
              "forcats", "openxlsx", "scales", "ggthemes", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Makes argparse object and arguments
parser <- ArgumentParser(description="Runs GSEA on gene-level genetic interaction scores.")
parser$add_argument("-i", "--input_file", type="character",
                    help="Path to file containing genetic interaction scores (rows) for set of queries (columns)")
parser$add_argument("-a", "--annotation_file", type="character", default=file.path("input", "Human_Reactome_October_01_2018_symbol.gmt"),
                    help="Path to gmt file containing pathway annotations [default %(default)s]")
parser$add_argument("-o", "--output_folder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("-q", "--query", type="character", default="all",
                    help="Query or queries (gene name or analysis ID; comma separated; case insensitive) to run GSEA [default %(default)s]")
parser$add_argument("--MIN_GENE", type="integer", default = 10,
                    help="Minimum number of genes to be considered in a pathway [default %(default)s]")
parser$add_argument("--MAX_GENE", type="integer", default = 500,
                    help="Maximum number of genes to be considered in a pathway [default %(default)s]")
parser$add_argument("--SIG_FDR", type="integer", default = 0.05,
                    help="FDR threshold to define significant pathway enrichment [default %(default)s]")
parser$add_argument("--SET_SEED", type="integer", default = 42,
                    help="Seed value to maintain result consistency [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

# Sets user parameters
score_file <- args$input_file
annotation_file <- args$annotation_file
output_folder <- args$output_folder
query <- strsplit(args$query, ",")[[1]]
MIN_GENE <- args$MIN_GENE
MAX_GENE <- args$MAX_GENE
SIG_FDR <- args$SIG_FDR
SET_SEED <- args$SET_SEED

# Check if required files exist
if (!file.exists(score_file)) {
  stop("Input gene-level score file does not exist")
}
if (!file.exists(annotation_file)) {
  stop("Pathway annotation file does not exist")
}

# Generate parent output folder if it doesn't already exist
if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Print user parameters to console
cat(sprintf("Hostname: %s\n", Sys.info()["nodename"]))
cat(sprintf("Start time: %s\n", format(Sys.time(), "%a %b %d %X %Y")))
cat(sprintf("Working directory: %s\nOutput directory: %s\n", getwd(), output_folder))
cat(sprintf("Gene-level genetic interaction score file: %s\n", score_file))
cat(sprintf("Pathway annotation file: %s\n", annotation_file))
cat(sprintf("GSEA parameters:\nPathway size: %i-%i\n\n", MIN_GENE, MAX_GENE))

######
# READ SCORES
######

# Read in score data
suppressWarnings({
  score <- fread(score_file, stringsAsFactors = FALSE, data.table = FALSE)
  rownames(score) <- score[,1]; score[,1] <- NULL
})
cat("---------------------------------------\n")
cat(sprintf("Read score data for %i genes across %i screens\n", nrow(score), ncol(score)))

# Subset data for queries of interest, if specified
if (length(which(query == "all")) == TRUE) {
  cat("Keeping all screens!\n")
  score_run <- score
} else {
  if (length(query) == 1) {
    cat(sprintf("Subsetting for screen %s\n", query))
    score_run <- score[,grep(query, colnames(score), ignore.case = TRUE), drop = FALSE]
  } else {
    cat(sprintf("* Subsetting for screens %s\n", paste(query, collapse=", ")))
    query_match <- paste(query, collapse="|")
    score_run <- score[,grep(query_match, colnames(score), ignore.case = TRUE), drop = FALSE]
  }
}

# Stop code if no data remains
if (!ncol(score_run)) {
  stop("No data left after subsetting (check spelling)\n")
} else {
  cat(sprintf("* Found %s screen(s) (%s)\n", ncol(score_run), paste(colnames(score_run), collapse=", ")))
}

######
# PATHWAY ANNOTATIONS
######

# Read in pathway file as a list
pathways <- gmtPathways(annotation_file)
pathways <- pathways[!duplicated(pathways)]

# Remove annotation source information from pathway names
names(pathways) <- gsub("\\%..*", "", names(pathways))

# Filter for pathway size (MIN_GENE, MAX_GENE)
pathway_size <- lapply(pathways, length)
pathways2 <- pathways[which(pathway_size > MIN_GENE & pathway_size < MAX_GENE)]

if (!length(pathways2)) {
  stop("No pathways left after size filtering! Choose new size filters")
} else {
  cat("\n---------------------------------------\n")
  cat(sprintf("Testing %s pathways from %s total (min size %s, max size %s)\n\n",
    length(pathways2), length(pathways), MIN_GENE, MAX_GENE))
}

# Get annotation source to add to file names
annotation_sources <- c("GO", "KEGG", "Reactome", "NetPath", "Corum")
annotation_name <- unlist(strsplit(basename(annotation_file), split = "_"))
annotation_source <- unique(grep(paste(annotation_sources, collapse = "|"),
    annotation_name, ignore.case = TRUE, value = TRUE))

if (!length(annotation_source)) {
  annotation_source <- "genericPathways"
}

######
# GSEA
######

# Loop GSEA run per screen
for (i in 1:ncol(score_run)) {

  # Define tested query
  query_i <- names(score_run)[i]
  cat(sprintf("\n[%i] Running GSEA on screen %s\n", i, query_i))

  # Prepare ranked data
  gi <- data.frame(gene = rownames(score_run), score = score_run[,i])
  gi <- setNames(gi$score, gi$gene) # vector format

  # Rank order data in decreasing order and remove missing values
  ## This will test the extreme tails of the ranked list of genes
  ## (ie strongest positive & negative interactions)
  gi_rank <- sort(gi, decreasing = TRUE)
  gi_rank <- gi_rank[!is.na(gi_rank)]

  # Run GSEA
  set.seed(SET_SEED)
  gsea <- fgsea(stats = gi_rank, pathways = pathways2, eps = 0)

  # Filtering for significant pathways
  gsea_sig <- as.data.table(filter(gsea, padj <= SIG_FDR))
  cat(sprintf("* Found %i enriched pathways (FDR <= %g)\n", nrow(gsea_sig), SIG_FDR))

  # Separate by negative / positive enrichment
  gsea_pos <- gsea_sig[NES > 0][order(NES, decreasing = TRUE),]
  gsea_neg <- gsea_sig[NES < 0][order(NES, decreasing = FALSE),]
  cat(sprintf("** positive enrichment: %s\n** negative enrichment: %s\n", nrow(gsea_pos), nrow(gsea_neg)))

  ######
  # VISUALIZATION
  ######

  if (!nrow(gsea_sig)) {
    cat("Nothing to see here...\n")
  } else {

    # Query-specific output folder to write results to
    query_folder <- paste(output_folder, query_i, sep="/")
    if (!dir.exists(query_folder)) { dir.create(query_folder) }

    # Main output file names
    out_file <- sprintf("%s/gsea_%s_%s", query_folder, query_i, annotation_source)

    ######
    # PATHWAY ENRICHMENT TABLES
    ######

    # Store results in list, if results exist
    # One sheet per set
    res <- list()
    if (nrow(gsea_pos)) {
      res[["positive"]] <- gsea_pos
    }
    if (nrow(gsea_neg)) {
      res[["negative"]] <- gsea_neg
    }

    # Write out results to file
    table_file <- sprintf("%s_pathway_table.xlsx", out_file)
    cat(sprintf("  => writing out enrichment results to %s.\n", basename(table_file)))
    write.xlsx(res, file = table_file)

    ######
    # PREPARE DATA FOR PLOTTING
    ######

    # Combine positive and negative results
    gsea_res <- c()
    if (nrow(gsea_pos)) {
      gsea_pos$sign <- "Positive"
      gsea_res <- rbind(gsea_res, gsea_pos)
    }
    if (nrow(gsea_neg)) {
      gsea_neg$sign <- "Negative"
      gsea_res <- rbind(gsea_res, gsea_neg)
    }

    # Get character vector of pathways
    pathway_res <- as.character(unlist(gsea_res[,"pathway"]))

    # Prepare data for heatmap
    gsea_heatmap <- gsea_res$leadingEdge
    names(gsea_heatmap) <- gsea_res$pathway

    # Melt to dataframe
    gsea_heatmap <- melt(gsea_heatmap)
    colnames(gsea_heatmap) <- c("Gene", "Pathway")
    gsea_heatmap$Gene <- as.character(gsea_heatmap$Gene)

    # Transform score vector to df and merge with gsea_heatmap
    gi_rank_df <- as.data.frame(gi_rank)
    gi_rank_df$Gene <- rownames(gi_rank_df)
    rownames(gi_rank_df) <- NULL
    colnames(gi_rank_df)[1] <- "gene_score"

    # Prevent merge from re-arranging columns using join
    gsea_heatmap <- left_join(gsea_heatmap, gi_rank_df, by = "Gene")

    # Add additional enrichment info
    gsea_cols <- gsea_res[,c("pathway", "padj", "NES", "size", "sign")]
    colnames(gsea_cols) <- c("Pathway", "FDR", "NES", "Size", "Sign")
    gsea_heatmap <- left_join(gsea_heatmap, gsea_cols, by = "Pathway")

    # Prevent ggplot from re-arranging pathway/gene levels
    gsea_heatmap$Pathway <- factor(gsea_heatmap$Pathway, levels = unique(gsea_res$pathway))
    gsea_heatmap$Gene <- factor(gsea_heatmap$Gene, levels = unique(gsea_heatmap$Gene))

    # Prepare data for dotplot
    # Grab top 20 positive / negative enrichments
    top_pos <- gsea_pos[order(gsea_pos$NES, decreasing = TRUE),][1:20,]
    top_neg <- gsea_neg[order(gsea_neg$NES, decreasing = FALSE),][1:20,]
    gsea_dotplot <- rbind(top_pos, top_neg)
    gsea_dotplot <- na.omit(gsea_dotplot)

    # Convert Sign to factor to prevent unwanted point re-arranging
    gsea_dotplot$sign <- factor(gsea_dotplot$sign, levels = unique(gsea_dotplot$sign))

    ######
    # LEADING EDGE PLOT
    ######

    # Draw out
    plot_file1 <- sprintf("%s_leadingEdge.pdf", out_file)
    cat(sprintf("  => plotting leading edge plot to %s.\n", basename(plot_file1)))

    pdf(plot_file1, width = 20, height = length(pathway_res)/2.5)
    plotGseaTable(
      pathways = pathways2[pathway_res],
      stats = gi_rank,
      fgseaRes = gsea_res,
      gseaParam = 0.5
    )
    invisible(dev.off())

    ######
    # HEATMAP
    ######

    # Draw out heatmaps
    plot_file2 <- sprintf("%s_heatmap.pdf", out_file)
    cat(sprintf("  => plotting enrichment heatmap to %s.\n", basename(plot_file2)))

    # Plot
    p2 <- ggplot(gsea_heatmap, aes(Gene, Pathway)) +
          facet_grid(. ~ Sign, scales = "free", space = "free") +
          geom_tile(aes(fill = gene_score), colour = "white") +
          labs(x = "Leading edge gene", y = "Pathway", title = query_i) +
          scale_fill_gradient2(low = muted("#386cb0"), high = "#ef3b2c") +
          theme_bw(base_size = 14) +
          theme(panel.background = element_rect(fill = "white"),
                panel.grid = element_blank(),
                plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
                axis.text.y = element_text(size = 8, hjust = 1, vjust = 1, face = "bold"),
                axis.title = element_text(size = 12),
                axis.ticks = element_blank(),
                legend.title = element_text(face = "bold", size = 10))

    ggsave(plot_file2, p2, width = length(unique(gsea_heatmap$Gene)) * 0.7,
           height = length(unique(gsea_heatmap$Pathway)) * 0.2, limitsize = FALSE)

    ######
    # DOTPLOT
    ######

    # Workaround to split long pathway names on multiple lines
    # (scales wrap_format() best suited for x axis labels)
    path <- gsea_dotplot$pathway
    path <- unlist(lapply(path, function(x) {
      x_split <- unlist(strsplit(x, ""))
      x_ins <- grep(" ", x_split)[which(grep(" ", x_split) > 50)][1]
      if (is.na(x_ins)) {
        return(x)
      } else {
        substr(x, x_ins, x_ins) <- "\n"
        return(x)
      }
    }))

    # Modify gsea_dotplot pathway names
    gsea_dotplot$pathway <- path

    # Set point fill
    p_fill <- unique(ifelse(gsea_dotplot$sign == "Negative", "cornflowerblue", "goldenrod1"))

    # Plot merged dot plot
    p3 <- gsea_dotplot %>%
      mutate(log10FDR = -log10(padj)) %>%
      ggplot(aes(x = log10FDR, y = fct_reorder(pathway, log10FDR), fill = sign)) +
        facet_grid("sign", scales = "free", space = "free") +
        geom_point(aes(size = abs(NES)), shape = 21, colour = "black") +
        #geom_vline(xintercept = -log10(SIG_FDR), linetype = "dotted", colour = "black", size = 0.75) +
        labs(y = NULL, x = "-log10(FDR)",
             title = sprintf("%s (%s)", query_i, annotation_source),
             fill = "Genetic interaction",
             size = "Normalized\nenrichment\nscore") +
        scale_fill_manual(values = p_fill, guide = FALSE) +
        theme_clean(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5),
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 10),
              legend.key.size = unit(0.1, "line"))

      # Estimate plot height
      #p_height <- length(unique(gsea_dotplot$pathway)) * 0.25

      # Draw out plot
      plot_file3 <- sprintf("%s_dotplot.pdf", out_file)
      cat(sprintf("  => plotting enrichment dotplot to %s.\n", basename(plot_file3)))
      ggsave(plot_file3, p3, width = 8, height = 9.5, limitsize = FALSE)
  }
}
