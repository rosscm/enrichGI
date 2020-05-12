# Running GO enrichment on gene-level pi scores using `enrichGO`

#' @param piScoreF (char) path to file with piScore data
#' @param fdrF (char) path to file with complementary FDR data
#' @param sigTailFDR (integer) set FDR threshold to select significant tails
#'		of Pi score distribution to test for GO enrichment (default=0.2)
#' @param keyType (char) keytype of input gene (default="SYMBOL")
#' @param OrgDb (char) genotype annotation data (default="org.Hs.eg.db")
#' @param ontGO (char) one of "MF", "BP", or "CC" sub-ontologies (default="BP")
#' @param pAdjMethod (char) one of "holm", "hochberg", "hommel", "bonferroni",
#'     "BH", "BY", "fdr", "none" (default="BH")
#' @param minGene (integer) min. number of genes to be considered in a pathway
#'    (default=10)
#' @param maxGene (integer) max. number of genes to be considered in a pathway
#'    (default=500)
#' @param pthres (integer) enrichGO pvalue cutoff (default=0.05)
#' @param qthres (integer) enrichGO qvalue cutoff (default=0.2)
#' @param nShowSig (integer) number of enriched pathways to display in output
#'    plots (default=10)

suppressMessages({
  require(enrichplot)
  require(clusterProfiler)
  require(ggplot2)
  require(cowplot)
  require(tidyr)
  require(xlsx)
})

source("../usefulFunctions.R")

# Paths to gene-level pi score data + pathway annotation files
dataDir <- "/Users/catherineross/GIN"
inDir   <- sprintf("%s/data/data_freeze_20190612", dataDir)

# Pi score data
piScoreF <- paste(inDir, grep("gi_PiScores_20.*.txt", list.files(inDir), value=T), sep="/")
fdrF     <- paste(inDir, grep("gi_FDR_20.*.txt", list.files(inDir), value=T), sep="/")

runEnrichGO <- function(piScoreF, fdrF,
                        sigTailFDR=0.2, keyType="SYMBOL",
                        OrgDb="org.Hs.eg.db",
                        ontGO="BP", pAdjMethod="BH",
                        minGene=10, maxGene=500,
                        pthres=0.05, qthres=0.2,
                        nShowSig=10) {

  dt <- format(Sys.Date(), "20%y%m%d")
  outDir <- sprintf("%s/res/%s_out_enrichGO_%s_sigTailFDR%i",
                    dataDir, dt, basename(inDir), sigTailFDR*100)

  if (file.exists(outDir)) {
    cat("\nOutput directory exists! Not overwriting\n\n")
    Sys.sleep(3)
  }

  if (!file.exists(outDir)) dir.create(outDir)

  ## WORK BEGINS
  sink(sprintf("%s/runEnrichGO.log", outDir))

  cat(sprintf("Hostname: %s\n", Sys.info()["nodename"]))
  cat(sprintf("Working directory: %s\nOutput directory: %s\n", getwd(), outDir))
  cat(sprintf("Start time: %s\n\n", format(Sys.time(), "%a %b %d %X %Y")))

  cat(sprintf("Gene-level Pi score file: %s\n", piScoreF))
  cat(sprintf("Gene-level FDR file: %s\n", fdrF))
  cat(sprintf("Pi score distribution tail threshold: %g\n\n", sigTailFDR))
  cat(sprintf("EnrichGO parameters:
  Key type = %s\n  Organism = %s\n  GO ontology = %s\n  p value cutoff = %g
  q value cutoff = %g\n  p value adjustment = %s\n  Pathway size = %i - %i\n",
  keyType, OrgDb, ontGO, pthres, qthres, pAdjMethod, minGene, maxGene))

  # Prepare data
  piScore <- read.delim(piScoreF, h=T, as.is=T)
  fdr     <- read.delim(fdrF, h=T, as.is=T)

  cat(sprintf("\nRead Pi score data for %i genes across %i query screens\n",
      nrow(piScore), ncol(piScore)))

  enrich_pos_all <- list()
  enrich_neg_all <- list()

  # Run enrichGO per screen
  for (i in seq_along(piScore)) {
    tryCatch({
      cat(sprintf("\n[%i] Running enrichGO on query screen %s...\n", i, names(piScore)[i]))

      # Read in GIs from pi score file (eg 17k library)
      dat_pi  <- setNames(piScore[,i], row.names(piScore))
      # Plus associated FDR values
      dat_fdr <- setNames(fdr[,i], row.names(fdr))

      # Filter GIs by specified FDR threshold (`sigTailFDR`)
      fdr_filter <- which(dat_fdr <= sigTailFDR)
      dat <- dat_pi[fdr_filter]

      # Get positive interactions
      pos_dat <- dat[which(dat > 0)]
      cat(sprintf("\n*Testing %g positive GIs (FDR <= %g) for GO enrichment",
                length(pos_dat), sigTailFDR))

      # Get negative interactions
      neg_dat <- dat[which(dat < 0)]
      cat(sprintf("\n*Testing %g negative GIs (FDR <= %g) for GO enrichment",
                length(neg_dat), sigTailFDR))

      # Run EnrichGO and generate gene-concept network plots
      cnetplot_list <- list(p1=NULL, p2=NULL)
      heatplot_list <- list(q1=NULL, q2=NULL)

      outF <- sprintf("%s/enrichGO_%s_tailFDR%g_%s",
                  outDir, basename(piScoreF), sigTailFDR*100, names(piScore)[i])

      plotRes <- function(df, resSign) {

        res <- enrichGO(gene=names(df),
                        universe=names(piScore),
                        keyType=keyType,
                        OrgDb=OrgDb,
                        ont=ontGO,
                        pAdjustMethod=pAdjMethod,
                        minGSSize=minGene,
                        maxGSSize=maxGene,
                        pvalueCutoff=pthres,
                        qvalueCutoff=qthres)

        # Write out results summary of all queries to object
        if (resSign == "pos") {
          if (nrow(res) == 0) {
            enrich_pos_all[[i]] <<- as.data.frame(NA)
          } else {
            enrich_pos_all[[i]] <<- as.data.frame(res$Description)
          }
          names(enrich_pos_all[[i]]) <<- names(piScore)[i]
        }
        if (resSign == "neg") {
          if (nrow(res) == 0) {
            enrich_neg_all[[i]] <<- as.data.frame(NA)
          } else {
            enrich_neg_all[[i]] <<- as.data.frame(res$Description)
          }
          names(enrich_neg_all[[i]]) <<- names(piScore)[i]
        }

        # Calls GOSemSim to calculate similarities among GO terms and remove those
        # highly similar terms by keeping one representative term
        #cat("Simplifying results via GoSemSim\n")
        #res_sim <- simplify(res, cutoff=0.7, by="p.adjust", select_fun=min)

        # Skip downstream steps if query has no significant enrichment
        if (nrow(res) == 0) {
          cat(sprintf("\nNo significant %s enrichment!\n\n", resSign))
          return(NULL)
        }

        cat(sprintf("\n\n**Got %g enriched GO terms (%s)", nrow(res), resSign))

        # Print out GO enrichment results
        cat(sprintf("\n  => writing out enrichment results to %s_%s.txt",
              basename(outF), resSign))
        write.table(as.data.frame(res),
                    file=sprintf("%s_%s.txt", outF, resSign),
                    col=T, row=F, quote=F, sep="\t")

        # Gene-concept network plots
        # NOTE suppressing ggplot messages about replacing colour + fill scales
        ## TRY cnetplot(fixed=F) to adjust network manually (2018-11-30)
        suppressMessages({
          p <- cnetplot(res,
                        categorySize="geneNum",
                        showCategory=nShowSig,
                        foldChange=dat) +
                  ggplot2::ggtitle(sprintf("GO enrichment (%s)", resSign)) +
                  ggplot2::scale_colour_gradient2(low="blue",
                                                  high="red",
                                                  midpoint=0,
                                                  name="Pi score")

          # Heat maps
          q <- heatplot(res,
                        showCategory=nShowSig,
                        foldChange=dat) +
                 ggplot2::ggtitle(sprintf("GO enrichment (%s)", resSign)) +
                 ggplot2::theme(plot.title=element_text(hjust=0.5)) +
                 ggplot2::scale_fill_gradient2(low="blue",
                                               high="red",
                                               midpoint=0,
                                               name="Pi score")
        })

        if (resSign == "pos") {
          cnetplot_list$p1 <<- p
          heatplot_list$q1 <<- q
        }
        if (resSign == "neg") {
          cnetplot_list$p2 <<- p
          heatplot_list$q2 <<- q
        }
     }

     plotRes(df=pos_dat, resSign="pos")
     plotRes(df=neg_dat, resSign="neg")

     # Write out plots
     # cnetplot
     plotF <- sprintf("%s_cnetplot", outF)
     cat(sprintf("\n  => plotting gene-concept network plot to %s_cnetplot.pdf",
          basename(outF)))

     pdf(sprintf("%s.pdf", plotF), width=25, height=12)
     p_all <- plot_grid(cnetplot_list$p1, cnetplot_list$p2)
     title <- ggdraw() + draw_label(names(piScore)[i], fontface="bold")
     print(plot_grid(title, p_all, ncol=1, rel_heights=c(0.1, 1)))
     invisible(dev.off())

     # heatplot
     plotF <- sprintf("%s_heatmap", outF)
     cat(sprintf("\n  => plotting heat map to %s_heatmap.pdf\n", basename(outF)))

     pdf(sprintf("%s.pdf", plotF), width=22, height=10)
     q_all <- plot_grid(heatplot_list$q1, heatplot_list$q2, ncol=1, nrow=2)
     title <- ggdraw() + draw_label(names(piScore)[i], fontface="bold")
     print(plot_grid(title, q_all, ncol=1, rel_heights=c(0.1, 1)))
     invisible(dev.off())

   }, error = function(e){cat("\n\nERROR:", conditionMessage(e), "\n")})
  }

   # Get enrichment summaries; excel sheets with all enriched pathways per query
   # Summary of enrichment for each query (n cols = n queries)
   pos_all <- do.call("cbind.na", enrich_pos_all)
   write.xlsx(as.data.frame(pos_all),
              file=sprintf("%s/%s_summary_pos_enrich_query.xlsx",
                           outDir, basename(inDir)),
              row.names=F, showNA=F)

   neg_all <- do.call("cbind.na", enrich_neg_all)
   write.xlsx(as.data.frame(neg_all),
              file=sprintf("%s/%s_summary_neg_enrich_query.xlsx",
                           outDir, basename(inDir)),
              row.names=F, showNA=F)

    # Pathway-level summary across all queries (number of times a pathway is
    # enriched across n queries)
    pos_all_col <- gather(pos_all)
    pos_tally <- as.data.frame(table(pos_all_col$value))
    pos_tally <- pos_tally[order(pos_tally$Freq, decreasing=T),]
    colnames(pos_tally)[1] <- "Pathway"
    write.xlsx(pos_tally,
               file=sprintf("%s/%s_summary_pos_enrich_tally.xlsx",
                            outDir, basename(inDir)),
               row.names=F, showNA=F)

    neg_all_col <- gather(neg_all)
    neg_tally <- as.data.frame(table(neg_all_col$value))
    neg_tally <- neg_tally[order(neg_tally$Freq, decreasing=T),]
    colnames(neg_tally)[1] <- "Pathway"
    write.xlsx(neg_tally,
               file=sprintf("%s/%s_summary_neg_enrich_tally.xlsx",
                            outDir, basename(inDir)),
                row.names=F, showNA=F)

  cat(sprintf("\nEnd time: %s\n\n", format(Sys.time(), "%a %b %d %X %Y")))
  sink()
}

runEnrichGO(piScoreF=piScoreF, fdrF=fdrF, sigTailFDR=0.5,
            keyType="SYMBOL", OrgDb="org.Hs.eg.db",
            ontGO="BP", pAdjMethod="BH",
            minGene=10, maxGene=500,
            pthres=0.05, qthres=0.05,
            nShowSig=15)
