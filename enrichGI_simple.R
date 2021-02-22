# Loads packages
packages <- c("fgsea", "data.table", "dplyr", "reshape2", "ggplot2", "scales",
              "forcats", "openxlsx", "Hmisc", "ggthemes", "argparse", "rio", "purrr")

for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# ============================= #
# User inputs (change these)
input_file <- "~/projects/GIN/data/forKuheli/CHEM022_032_drugz.xlsx"
annotation_file <- "~/projects/GIN/anno/Human_GO_bp_no_GO_iea_October_01_2018_symbol.gmt"
output_folder <- sprintf("~/projects/GIN/res/%s_out_forKuheli_GSEA", format(Sys.Date(), "20%y%m%d"))

# Parameters
gene_col = "gene" # column with gene symbols
score_col = "normZ" # column with gene scores
MIN_GENE = 10 # min number of genes in pathway
MAX_GENE = 300 # max number of genes in pathway
SIG_FDR = 0.05 # enrichment FDR significance
# ============================= #

# Read input data
input <- lapply(input_file, import_list)
input <- flatten(input)

# Generate output directory
if (!file.exists(output_folder)) dir.create(output_folder)

# Read in pathway file as a list and clean pathway names
pathways <- gmtPathways(annotation_file)
names(pathways) <- gsub("\\%..*", "", names(pathways))
pathway_size <- lapply(pathways, length)
pathways2 <- pathways[which(pathway_size >= MIN_GENE & pathway_size <= MAX_GENE)]

# Get annotation source to add to file names
annotation_sources <- c("GO", "KEGG", "Reactome", "NetPath", "Corum")
annotation_name <- unlist(strsplit(basename(annotation_file), split = "_"))
annotation_source <- unique(grep(paste(annotation_sources, collapse = "|"),
    annotation_name, ignore.case = TRUE, value = TRUE))

# Loop GSEA run per screen
for (i in seq_along(input)) {

  # Define tested screen
  screen <- input[[i]]
  screen_i <- names(input)[i]
  cat(sprintf("\n[%i] Running GSEA on screen %s\n", i, screen_i))

  # Prepare ranked data
  gi <- data.frame(gene = screen[,gene_col], score = screen[,score_col])
  gi <- setNames(gi$score, gi$gene) # vector format

  # Rank order data in decreasing order and remove missing values
  ## This will test the extreme tails of the ranked list of genes
  ## (ie strongest positive & negative interactions)
  gi_rank <- sort(gi, decreasing = TRUE)
  gi_rank <- gi_rank[!is.na(gi_rank)]

  # Run GSEA
  set.seed(42)
  gsea <- fgsea(stats = gi_rank, pathways = pathways2)

  # Filtering for significant pathways
  gsea_sig <- as.data.table(filter(gsea, padj <= SIG_FDR))
  n_sig <- nrow(gsea_sig)
  cat(sprintf("* Found %i enriched pathways (FDR <= %g)\n", n_sig, SIG_FDR))

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

    # Screen-specific output folder to write results to
    screen_folder <- paste(output_folder, screen_i, sep="/")
    if (!dir.exists(screen_folder)) { dir.create(screen_folder) }

    # Main output file names
    out_file <- sprintf("%s/gsea_%s_%s", screen_folder, screen_i, annotation_source)

    ######
    # PATHWAY ENRICHMENT TABLES
    ######

    # Store results in list, if results exist
    # One sheet per set
    res <- list()
    if (nrow(gsea_pos)) { res[["positive"]] <- gsea_pos }
    if (nrow(gsea_neg)) { res[["negative"]] <- gsea_neg }

    # Write out results to file
    table_file <- sprintf("%s_pathway_table.xlsx", out_file)
    cat(sprintf("  => writing out enrichment results to %s.\n", basename(table_file)))
    write.xlsx(res, file = table_file)

    ######
    # LEADING EDGE PLOT
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

    # Add column to gsea_res of number of leadingEdge genes in pathway
    gsea_res$count <- ""
    for (i in 1:nrow(gsea_res)) {
      gsea_res$count[i] <- length(as.character(unlist(gsea_res[i,"leadingEdge"])))
    }

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
    # DOTPLOT
    ######

    # Prepare data
    gsea_plot <- gsea_res$leadingEdge
    names(gsea_plot) <- gsea_res$pathway

    # Melt to dataframe for ggplot
    gsea_plot <- reshape2::melt(gsea_plot)
    colnames(gsea_plot) <- c("Gene", "Pathway")
    gsea_plot$Gene <- as.character(gsea_plot$Gene)

    # Transform score vector to df and merge with gsea result df
    gi_rank_df <- as.data.frame(gi_rank)
    gi_rank_df$Gene <- rownames(gi_rank_df)
    rownames(gi_rank_df) <- NULL
    colnames(gi_rank_df)[1] <- "gene_score"

    # Prevent merge from re-arranging columns using join (2018-11-08)
    final <- left_join(gsea_plot, gi_rank_df, by = "Gene")

    # Add additional enrichment info
    gsea_sign <- gsea_res[,c("pathway", "padj", "NES", "size", "count", "sign")]
    colnames(gsea_sign) <- c("Pathway", "FDR", "NES", "Size", "Count", "Sign")
    final <- left_join(final, gsea_sign, by = "Pathway")

    # Prevent ggplot from re-arranging pathway/gene levels
    final$Pathway <- factor(final$Pathway, levels = unique(gsea_res$pathway))
    final$Gene <- factor(final$Gene, levels = unique(final$Gene))

    # Remove Gene and score_score columns from final to get unique set of pathways
    final2 <- unique(final[,-c(1,3)])
    final2$Pathway <- as.character(final2$Pathway)

    # Calculate geneRatio
    final2$geneRatio <- as.numeric(final2$Count) / as.numeric(final2$Size)

    # Uncapitalize pathway names
    final2$Pathway <- tolower(final2$Pathway)
    final2$Pathway <- capitalize(final2$Pathway)

    # Convert Sign to factor format to prevent unwanted point re-arranging
    final2$Sign <- factor(final2$Sign, levels = unique(final2$Sign))

    # Grab top 20 positive / negative enrichments
    top_pos <- final2[order(final2$NES, decreasing = TRUE),][1:20,]
    top_neg <- final2[order(final2$NES, decreasing = FALSE),][1:20,]
    top_plot <- rbind(top_pos, top_neg)

    # Set point fill
    p_fill <- unique(ifelse(top_plot$Sign == "Negative", "cornflowerblue", "goldenrod1"))

    # Set x-axis limits so points aren't cut off from plot window
    xmin <- floor(min(top_plot$NES))
    xmax <- ceiling(max(top_plot$NES))

    # Plot merged dot plot
    p <- ggplot(top_plot, aes(x = NES,
                              y = fct_reorder(Pathway, NES),
                              fill = Sign)) +
            facet_grid("Sign", scales = "free", space = "free") +
            geom_point(aes(size = -log10(FDR)), shape = 21, colour = "black") +
            geom_vline(xintercept = 0, linetype = "dotted", colour = "black", size = 0.75) +
            labs(y = NULL, x = "Normalized enrichment score",
                 title = sprintf("%s (%s)", screen_i, annotation_source),
                 fill = "Genetic interaction",
                 size = "-log10(FDR)") +
            xlim(xmin, xmax) +
            scale_fill_manual(values = p_fill) +
            theme_few(base_size = 14) +
            theme(plot.title = element_text(hjust = 0.5),
                  legend.title = element_text(size = 10),
                  legend.text = element_text(size = 10),
                  legend.key.size = unit(0.1, "line"),
                  strip.text.y = element_blank(),
                  panel.grid.major.y = element_line(linetype = "dotted", colour = "lightgrey"))

      # Draw out plot
      plot_file <- sprintf("%s_dotplot_merged.pdf", out_file)
      cat(sprintf("  => plotting enrichment dotplot to %s.\n", basename(plot_file)))
      ggsave(plot_file, p, width = 10, height = 12.5, limitsize = FALSE)
  }
}
