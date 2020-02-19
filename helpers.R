######################################################################################
################################## Helper Functions ##################################
######################################################################################

# Pretty timestamp
timestamp2 <- function() {
  # Partially from https://stackoverflow.com/questions/1962278/dealing-with-timestamps-in-r
  now <- Sys.time()
  timeList <- unclass(as.POSIXlt(now))
  secNow <- ifelse(round(timeList$sec) < 10, paste0(0, round(timeList$sec)), round(timeList$sec))
  paste0("[", timeList$hour,":",timeList$min,":",secNow,
      " ", (timeList$mon+1), "/", timeList$mday, "/", (timeList$year + 1900), "]", sep = "")
}

# Convenient wrapper for Trendy and GSEA
doTrendyAnalysis <- function(countsNorm, timeVec, cores = 45) {
  res <- trendy(Data = countsNorm, NCores = cores,
                tVectIn = timeVec, maxK = 2)
  resNow <- results(res)
  res.top <- topTrendy(resNow, adjR2Cut = 0)
  ranks <- res.top$AdjustedR2 * res.top$Segment.Trends[,c(1)]
  names(ranks) <- names(res.top$AdjustedR2)
  ranks <- ranks[order(ranks, decreasing = F)]
  pltDF <- data.frame(
    "geneName" = names(ranks),
    "value" = ranks,
    "rank" = seq(1, length(ranks)),
    stringsAsFactors = F
  )
  TERM2GENE <- correlationAnalyzeR::getTERM2GENE("complex")
  GSEARes <- correlationAnalyzeR::myGSEA(ranks, TERM2GENE)
  GSEARes[["rankDF"]] <- pltDF
  return(GSEARes)
}

#Simple function for ranging a vector from 0 to 1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Categorizes Bulk RNA-Seq metadata using custom regex dictionary
categorizeMetaData <- function(metadata, cols, dictionary) {
  for (i in 1:length(names(dictionary))) {
    termNow <- names(dictionary)[i]
    dictNow <- dictionary[[i]]
    posInd <- c()
    negInd <- c()
    metadata$newColNow <- 0
    if ("yes" %in% names(dictNow)) {
      for (j in 1:length(cols)) {
        colNow <- cols[j]
        colIndNow <- which(colnames(metadata) == colNow)
        posInd <- c(posInd, grep(x = metadata[,colIndNow], 
                                 pattern = paste0(dictNow$yes, collapse = "|"), 
                                 perl = T, ignore.case = T))
        if (is.null(dictNow$no)) {
          negInd <- c(negInd)
        } else {
          negInd <- c(negInd, grep(x = metadata[,colIndNow], pattern = paste0(dictNow$no, collapse = "|"), perl = T, ignore.case = T))
        }
      }
      posInd <- unique(posInd[! posInd %in% negInd])
    } else {
      for (j in 1:length(cols)) {
        colNow <- cols[j]
        colIndNow <- which(colnames(metadata) == colNow)
        posInd <- c(posInd, grep(x = metadata[,colIndNow], pattern = paste0(dictNow, collapse = "|"), perl = T, ignore.case = T))
      }
    }
    metadata$newColNow[posInd] <- 1
    colnames(metadata)[which(colnames(metadata) == "newColNow")] <-  termNow
  }
  return(metadata)
}

# Finds number of samples gene is expressed in
nonZeroSamps <- function(row) {
  return(sum(row > 0))
}

# Convert to ScanPy Format
makeH5AD <- function(srt) {
  if (! dir.exists("Data/scRNASeq/h5ADs/")) {
    dir.create("Data/scRNASeq/h5ADs/")
  }
  scanpy <- reticulate::import("scanpy")
  datNow <- Matrix::t(srt@assays$RNA@counts)
  studyNow <- unique(srt$SRA)
  obsNow <- data.frame(
    Study = srt$SRA, 
    Sample = srt$SRS,
    Condition = srt$TissueSite,
    Chemistry = srt$chemistry,
    Counter = srt$Counter,
    Species = srt$Species,
    row.names = c(0:(length(srt$orig.ident) -1))
  )
  dir.create(paste0("Data/scRNASeq/h5ADs/", studyNow))
  annData <- scanpy$AnnData(X = datNow)
  annData$write_h5ad(filename = paste0("Data/scRNASeq/h5ADs/", studyNow, "/", studyNow, ".h5ad"))
  write.table(obsNow, sep = "\t", quote = F, row.names = T,
              file = paste0("Data/scRNASeq/h5ADs/", studyNow, "/", studyNow, "_cellMetaInfo.txt"))
  genesNow <- data.frame(
    geneName = rownames(srt), row.names = c(0:(length(rownames(srt)) -1))
  )
  write.table(genesNow, sep = "\t", quote = F, row.names = T,
              file = paste0("Data/scRNASeq/h5ADs/", studyNow, "/", studyNow, "_genesNames.txt"))
}

# Automatically downsample a Seurat object
downsampleSRT <- function(srt, n = 5000, sampleCol = "SRA") {
  colID <- which(colnames(srt@meta.data) == sampleCol)
  keepIDs <- unique(unlist(lapply(unique(srt@meta.data[,colID]), FUN = function(x) {
    toSamp <- names(srt$orig.ident)[srt@meta.data[,colID] == x]
    if (length(toSamp) < n) {
      toSamp
    } else {
      sample(toSamp, size = n)
    }
  })))
  srtNow <- subset(srt, cells = keepIDs)
  return(srtNow)
}

# Apply relative and absolute QC filtering cutoffs to Seurat object
autoQCSRT <- function(srt, species = "human", whichColumn = "SRA",
                      maxPercentMT = 20, 
                      minNumGenes = 200, minNumCounts = 1000,
                      maxNumGenes = Inf, maxNumCounts = Inf) {
  
  srtMeta <- srt@meta.data
  colId <- which(colnames(srtMeta) == whichColumn)
  origCount <- length(srtMeta[,colId])
  if (! "percent.mt" %in% colnames(srtMeta) |
      any(is.na(srtMeta$percent.mt))) {
    cat("\nCalculating percent MT\n")
    if (species == "human") {
      pattern <- "^MT-"
    } else {
      pattern <- "^mt-"      
    }
    srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = pattern)
    srtMeta <- srt@meta.data
  }
  # Filter srt based on absolute thresholds
  print(length(srt$orig.ident))
  cellsKeep <- rownames(srtMeta)[which(
    srtMeta$nCount_RNA >= minNumCounts &
      srtMeta$nCount_RNA <= maxNumCounts &
      srtMeta$nFeature_RNA >= minNumGenes &
      srtMeta$nFeature_RNA <= maxNumGenes &
      srtMeta$percent.mt <= maxPercentMT
  )]
  srt <- subset(srt, cells = cellsKeep )
  srtMeta <- srt@meta.data
  outCellsList <- lapply(unique(srtMeta[,colId]), FUN = function(x) {
    metaNow <- srtMeta[srtMeta[,colId] == x,]
    # Remove statistical outliers
    outlier1 <- boxplot.stats(metaNow$nCount_RNA)$out
    outlier2 <- boxplot.stats(metaNow$nFeature_RNA)$out
    outlier3 <- boxplot.stats(metaNow$percent.mt)$out
    outCells <- rownames(metaNow)[metaNow$nCount_RNA %in% outlier1 |
                                    metaNow$nFeature_RNA %in% outlier2 |
                                    metaNow$percent.mt %in% outlier3]
    outCells
  })
  badCellList <- unique(unlist(outCellsList, use.names = FALSE))
  goodCells <- rownames(srtMeta)[! rownames(srtMeta) %in% badCellList]
  srt <- subset(srt, cells = goodCells)
  newCount <- length(srtMeta[,colId])
  percentRemaining <- round(100*(newCount/(origCount)))
  cat(percentRemaining, "percent of cells passed QC!\n")
  return(srt)
}

# Convenience function to make feature plots similar to Seurat
featurePlot2 <- function(pltData, feature, 
                         featureType = c("column", "gene", "pathway"),
                         reduction = "umap") {
  require(ggpubr)
  
  
  
  
} 

# Makes GSEA plots in the paper style
gseaPlot2 <- function(EGMT, ID, title, scoreName = "Ranking metric") {
  # title = "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
  # ID = "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS"
  # EGMT <- GSEAResEWS$EGMT
  # scoreName = "Ranking metric"
  
  gp1 <- enrichplot::gseaplot2(EGMT, title = title,
                               subplots = c(1), 
                               geneSetID = ID)
  maxES <- gp1$data$x[which.max(abs(gp1$data$runningScore))]
  if (gp1$data$runningScore[which.max(abs(gp1$data$runningScore))] > 0) {
    limits <- c(0, (max(gp1$data$runningScore)+
            (max(gp1$data$runningScore)*.05)))
  } else {
    limits <- c((min(gp1$data$runningScore)+
                      (min(gp1$data$runningScore)*.05)), 0)
  }
  gp1 <- gp1 + 
    theme_pubr(border = T, base_size = 22) + 
    geom_vline(xintercept = maxES, colour = "red", linetype = 2) +
    ylab("Enrichment score") + theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    rremove("legend") + rremove("x.ticks") + rremove("x.text") +
    scale_y_continuous(expand = c(0, 0), limits = limits)
  
  gp2 <- enrichplot::gseaplot2(EGMT, title = NULL,
                               subplots = c(3), 
                               geneSetID = ID)
  gp2 <- gp2 + 
    theme_pubr(border = T, base_size = 22) + 
    ylab(scoreName) + theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    rremove("legend") + rremove("x.ticks") + rremove("x.text")
  ga <- ggarrange(gp1, gp2, ncol = 1, nrow = 2)
  ga
}

# Convenient wrapper for DimPlots in this publication
DimPlot2 <- function(data, mapping, plotName, height = 7.5, width = 8,
                     ptSize = .4,
                     legend = TRUE, silent = TRUE) {
  g1 <- ggplot(data = data, mapping = mapping) + geom_point(size = ptSize) +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
  ggsave(g1, filename = paste0("Figures/", plotName, ".png"), height = height, 
         width = width)
  if (legend) {
    g1 <- ggplot(data = data, mapping = mapping) + geom_point(size = ptSize) + theme_pubr(border = T, base_size = 22) +
      guides(colour = guide_legend(override.aes = list(size=3))) + 
      theme(legend.position="right") +
      theme(#border = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
    ggsave(g1, filename = paste0("Figures/", plotName, "_legend.png"), height = height, 
           width = round((width * 1.33), digits = 2))
  }
  if (! silent) {
    return(g1)
  } 
}

# Convenience function for PHATE correlation plots
doPhateCorrPlot <- function(nameNow, resSamples) {
  phateRes <- resSamples[[nameNow]]
  pltDF <- phateRes$rankDF
  pltDF$gene <- ""
  m <- length(pltDF$rank)-4
  n <- max(pltDF$value)
  o <- min(pltDF$value)
  pltDF$gene[(pltDF$geneName %in% c("FLI1", "FANCI", "BRCA1", "FEN1") |
                pltDF$rank < 4 | pltDF$rank > m)] <-as.character(pltDF$geneName)[(pltDF$geneName %in% c("FLI1", "FANCI", "BRCA1", "FEN1")|
                                                                                    pltDF$rank < 4 | pltDF$rank > m)] 
  
  
  pltDFLab <- pltDF[pltDF$gene != "",]
  g1 <- ggplot(pltDF, mapping = aes_string(x = "rank", y = "value")) +
    geom_point(size = .8, color = ifelse(pltDF$gene == "", "black", "red")) + 
    
    geom_hline(yintercept = 0, linetype="dashed", color = "grey") +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    theme_pubr(border = T, base_size = 22) +
    scale_x_continuous(limits=c((m*-1.05 + m), (m*1.3))) +
    scale_y_continuous(limits=c((o*1.1), (n*1.1))) +
    geom_text_repel( mapping = aes_string(x = "rank", 
                                          y = "value",
                                          label = "gene"), min.segment.length = .1,
                     data = pltDFLab, point.padding = 1, seed = 42, nudge_y = (o*.075), nudge_x = 3000) +
    ylab(paste0(nameNow, " correlation score")) +
    rremove("xlab") +
    rremove("legend")
  return(g1)
}

# Convenience wrapper for msigdbr
getTERM2GENE <- function(GSEA_Type = c("simple"),
                         Species = c("hsapiens", "mmusculus"),
                         sampler = FALSE) {
  
  # Species = "hsapiens"
  # GSEA_Type = "simple"
  # sampler = FALSE
  
  if (Species == "hsapiens") {
    msigSpec <- "Homo sapiens"
  } else {
    msigSpec <- "Mus musculus"
  }
  
  # Get data object
  MDF <- msigdbr::msigdbr(species = msigSpec)
  MDF$gs_subcat <- gsub(MDF$gs_subcat, pattern = "CP:", replacement = "", perl = TRUE)
  MDF$gs_cat <- paste0(MDF$gs_cat, ":", MDF$gs_subcat)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = ":$", replacement = "", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C1", replacement = "Cytogenic bands", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C6", replacement = "Oncogenic signatures", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C7", replacement = "Immunological signatures", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C2:", replacement = "", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C5", replacement = "GO", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "H", replacement = "Hallmark", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CP", replacement = "Canonical pathways", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "CGP", replacement = "Perturbations", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CGN", replacement = "Cancer gene neighborhoods", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C4:CM", replacement = "Cancer modules", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:MIR", replacement = "miRNA targets", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "C3:TFT", replacement = "TF targets", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "BIOCARTA", replacement = "BioCarta", perl = TRUE)
  MDF$gs_cat <- gsub(MDF$gs_cat, pattern = "REACTOME", replacement = "Reactome", perl = TRUE)
  
  # Filter for pathways of interest
  optionsNow <- c("simple", "complex", unique(MDF$gs_cat))
  if (! all(GSEA_Type %in% optionsNow)) {
    stop("\nPlease enter a valid GSEA_Type. Use ?getTERM2GENE to see available options.\n")
  }
  categories <- c()
  if ("simple" %in% GSEA_Type) {
    categories <- c(categories, "Hallmark", "Perturbations", "BioCarta",
                    "GO:BP", "KEGG", "Canonical pathways", "Reactome", "GO:MF", "GO:CC", "PID")
  }
  if ("complex" %in% GSEA_Type) {
    categories <- c(categories, optionsNow)
  }
  categories <- unique(c(categories, GSEA_Type))
  TERM2GENE <- MDF %>%
    filter(.data$gs_cat %in% categories) %>%
    select(.data$gs_name, .data$gene_symbol)
  
  
  if (sampler) {
    print("Using sampler!")
    set.seed(1)
    TERM2GENE <- TERM2GENE[sample(nrow(TERM2GENE), size = 100000),]
  }
  return(TERM2GENE)
}

# Convenience wrapper for making supplemental cluster marker plots
doClusterMarkerPlot <- function(srt, srtMarkers, outName) {
  # # Bug testing
  # srt <- BDInt
  # srtMarkers <- BDIntMarkers
  # outName = "Figures/Fig1Cxiii"
  
  dir.create(paste0("Figures/", outName), showWarnings = F)
  srtMarkers <- srtMarkers[srtMarkers$p_val_adj < .05,]
  srtMarkers$p_val_adj[srtMarkers$p_val_adj == 0] <- .Machine$double.xmin
  clusterList <- unique(as.numeric(srt$seurat_clusters))
  clusterList <- c(0, clusterList[order(clusterList)])
  names(clusterList) <- paste0("Cluster: ", clusterList)
  pltDat <- srt@meta.data
  pltDat$UMAP_1 <- srt@reductions$umap@cell.embeddings[,c(1)]
  pltDat$UMAP_2 <- srt@reductions$umap@cell.embeddings[,c(2)]
  TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple")
  for (i in 1:length(clusterList)) {
    clusterName <- names(clusterList)[i]
    print(clusterName)
    clusterNow <- clusterList[clusterName]
    pltDat$Selected <- FALSE
    pltDat$Selected[pltDat$seurat_clusters == clusterNow] <- TRUE
    pltDat$Selected <- factor( pltDat$Selected, levels = c(TRUE, FALSE))
    pltDat <- pltDat[order(pltDat$Selected, decreasing = TRUE),]
    g1 <- ggplot(data = pltDat, mapping = aes_string(x = "UMAP_1",
                                                     y = "UMAP_2",
                                                     color = "Selected")) + 
      geom_point(size = .6) + theme_pubr(border = T, base_size = 22) +
      guides(colour = guide_legend(override.aes = list(size=3))) + 
      theme(legend.position="right") +
      labs(title = clusterName) +
      theme(#border = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
      scale_color_manual(values=c("#CC0909", "#C4BFBF")) +
      rremove("legend")
    wordData <- srtMarkers[srtMarkers$cluster == clusterNow,]
    wordData <- data.frame(
      "word" = wordData$gene,
      "freq" = -log10(wordData$p_val_adj) * sign(wordData$avg_logFC),
      stringsAsFactors = FALSE
    )
    posCloud <- wordData[wordData$freq > 0,]
    posCloud <- posCloud %>% top_n(n = 6, wt = freq)
    negCloud <- wordData[wordData$freq < 0,]
    negCloud$freq <- -1*(negCloud$freq)
    negCloud <- negCloud %>% top_n(n = 6, wt = freq)
    if (! length(posCloud$word) | ! length(negCloud$word)) {
      warning("Not enough markers for cluster ", clusterNow)
      next
    }
    pal <-  colorRampPalette(brewer.pal(9,"Reds"))(length(posCloud$word))
    plot.new()
    dev.off()
    posCloud <- posCloud[order(posCloud$freq, decreasing = T),]
    posCloud$freq <- order(posCloud$freq, decreasing = F)
    wordcloud(words = posCloud$word, max.words = 6,
              random.order = F,  scale=c(2.5,.1),
              fixed.asp = T, rot.per = 0,
              freq = posCloud$freq, colors = pal)
    gridGraphics::grid.echo()
    pt <- grid.grab()
    g2 <- as.ggplot(pt)
    dev.off()
    pal <-  colorRampPalette(brewer.pal(9,"Greys"))(length(negCloud$word))
    plot.new()
    dev.off()
    negCloud <- negCloud[order(negCloud$freq, decreasing = T),]
    negCloud$freq <- order(negCloud$freq, decreasing = F)
    wordcloud(words = negCloud$word, max.words = 6,
              random.order = F, scale=c(2.5,.1),
              rot.per = 0, fixed.asp = T,
              freq = negCloud$freq, colors = pal)
    gridGraphics::grid.echo()
    pt <- grid.grab()
    g3 <- as.ggplot(pt)
    dev.off()
    
    # Pathway plots
    posGenes <- wordData$word[wordData$freq > 0]
    negGenes <- wordData$word[wordData$freq < 0]
    g4 <- doPathEnrichPlot(posGenes, negGenes, 
                           TERM2GENE = TERM2GENE)
    
    # Arrange
    ggNow <- ggdraw() +
      draw_plot(g1, x = 0, y = .4, width = .4, height = .6) +
      draw_plot(g2, x = 0, y = 0, width = .2, height = .35) +
      draw_plot(g3, x = 0.25, y = 0, width = .2, height = .35) +
      draw_plot(g4, x = 0.47, y = 0, width = .53, height = 1) 
    ggsave(ggNow, file = paste0("Figures/", outName, "/Cluster_", clusterNow, ".png"),
           height = 9, width = 16)
  }
}

# Makes gene set IDs presentable for publications
fixStrings <- function(StringVec) {
  
  # # Bug testing
  # StringVec <- c("HALLMARK_APOPTOSIS", "GO_MIR21_TARGETS",
  #                "GSE121239_THING_HAPPENED", "CTTGAT_MIR381",
  #                "GTGCAGAG_EZH2", "GSE12309_WHATEVER",
  #                "MEISSNER_NPC_HCP_WITH_H3K4ME2",
  #                "BASSO_CD40_SIGNALING_DN",
  #                "GSE36078_UNTREATED_VS_AD5_INF_IL1R_KO_MOUSE_LUNG_DC_DN",
  #                "KEGG_ACUTE_MYELOID_LEUKEMIA",
  #                "GSE31082_CD4_VS_CD8_SP_THYMOCYTE_UP",
  #                "GSE3920_IFNA_VS_IFNG_TREATED_ENDOTHELIAL_CELL_UP",
  #                "GSE37605_FOXP3_FUSION_GFP_VS_IRES_GFP_TREG_C57BL6_UP",
  #                "GSE41176_UNSTIM_VS_ANTI_IGM_STIM_BCELL_24H_DN",
  #                "MORI_LARGE_PRE_BII_LYMPHOCYTE_UP",
  #                "RYAAAKNNNNNNTTGW_UNKNOWN",
  #                "GGGTGGRR_PAX4_03",
  #                "GGGNNTTTCC_NFKB_Q6_01",
  #                "AAAYWAACM_HFH4_01",
  #                "KANG_DOXORUBICIN_RESISTANCE_UP")
  
  StringVec <- gsub(StringVec, pattern = "_", replacement = " ")
  StringVec <- tolower(StringVec)
  StringVec <- stringr::str_to_title(StringVec)
  StringVec <- gsub(StringVec, pattern = "Iii", replacement = "III", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ii", replacement = "II", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Of ", replacement = " of ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " To ", replacement = " to ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " In ", replacement = " in ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " At ", replacement = " at ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " With ", replacement = " with ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Without ", replacement = " without ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Upon ", replacement = " upon ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " An ", replacement = " an ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " By ", replacement = " by ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " For ", replacement = " for ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Via ", replacement = " via ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lof ", replacement = " LOF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lof$", replacement = " LOF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Loh ", replacement = " LOH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Loh$", replacement = " LOH", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Arms ", replacement = " ARMS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Erms ", replacement = " ERMS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nadh ", replacement = " NADH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nadph ", replacement = " NADPH ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cll ", replacement = " CLL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cml ", replacement = " CML ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Aml ", replacement = " AML ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " All ", replacement = " ALL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Kim ALL", replacement = "Kim all", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "CMV ALL", replacement = "CMV all", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nhek ", replacement = " NHEK ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ner ", replacement = " NER ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nmda ", replacement = " NMDA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dc ", replacement = " DC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cd4 ", replacement = " CD4 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Cd8 ", replacement = " CD8 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gc ", replacement = " GC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hdl ", replacement = " HDL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dn$", replacement = " Down", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ldl ", replacement = " LDL ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Tcr ", replacement = " TCR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mdc ", replacement = " MDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bcr ", replacement = " BCR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Icp ", replacement = " ICP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hbv ", replacement = " HBV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Dlbcl", replacement = "DLBCL", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gist ", replacement = " GIST ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gist$", replacement = " GIST", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dp ", replacement = " DP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dn ", replacement = " DN ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " H2o2 ", replacement = " H2O2 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " With ", replacement = " with ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ntreg", replacement = "nTreg", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mlr ", replacement = " MLR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Gfp", replacement = "GFP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Vs ", replacement = " vs ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " And ", replacement = " and ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Wt ", replacement = " WT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ros$", replacement = " ROS", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ros ", replacement = " ROS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ko ", replacement = " KO ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Pdc ", replacement = " PDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Pdgf", replacement = "PDGF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Rna ", replacement = " RNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mrna", replacement = "mRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mirna", replacement = "miRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Sirna", replacement = "siRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Trna", replacement = "tRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ncrna", replacement = "ncRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Snrna", replacement = "snRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Rrna", replacement = "rRNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna$", replacement = "RNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "rna ", replacement = "RNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Flii", replacement = "Fli1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hcp ", replacement = " HCP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Tnf ", replacement = " TNF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Srp ", replacement = " SRP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Utr ", replacement = " UTR ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Dna ", replacement = " DNA ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Rdna", replacement = "rDNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hiv ", replacement = " HIV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hiv1 ", replacement = " HIV1 ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "dna", replacement = "DNA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Lps ", replacement = " LPS ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gmcsf ", replacement = " GMCSF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Gm Csf ", replacement = " GMCSF ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Balbc", replacement = "BALBc", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Lcmv", replacement = "LCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mcmv", replacement = "MCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Pcc ", replacement = " PCC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ecm", replacement = "ECM", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "G1s", replacement = "G1S", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G1 S ", replacement = " G1S ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G1 S", replacement = " G1S", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "G2m", replacement = "G2M", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G2 M ", replacement = " G2M ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " G2 M", replacement = " G2M", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Hcmv", replacement = "HCMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Pbmc", replacement = "PBMC", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Atp ", replacement = " ATP ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Atp$", replacement = " ATP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Gtp", replacement = "GTP", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Mut ", replacement = " MUT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Et Al", replacement = "et al", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cpg", replacement = "CPG", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Nkt ", replacement = " NKT ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Hsc ", replacement = " HSC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Ln ", replacement = " LN ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cmv", replacement = "CMV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bm ", replacement = " BM ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Bmdc ", replacement = " BMDC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Esc ", replacement = "ESC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Esc ", replacement = " ESC ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mcf10a", replacement = "MCF10A", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tca", replacement = "TCA", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Nkcell", replacement = "NK-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tcell", replacement = "T-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "T Cell", replacement = "T-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "B Cell", replacement = "B-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Hela", replacement = "HeLa", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Treg", replacement = "T-Reg", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Tconv", replacement = "T-Conv", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Bcell", replacement = "B-Cell", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Uv ", replacement = " UV ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = " Uv$", replacement = " UV", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gse", replacement = "GSE", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gnf2", replacement = "GNF2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Gcm", replacement = "GCM", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Morf", replacement = "MORF", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Kegg", replacement = "KEGG", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Pid", replacement = "PID", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^Go ", replacement = "GO ", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ Q[0-9]$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ Q[0-9] [0-9]+$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( [A-Za-z0-9]+ [0-9]+$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "^([GCATSMWNRYK][gcatsmwrnyk]+)( Unknown$)",
                    replacement = "\\U\\1\\E\\2", perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Mir([0-9]*.*)", replacement = "miR\\1",
                    perl = TRUE, ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Ifn([a-z])", perl = TRUE,
                    replacement = "IFN\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Il([0-9]+)", perl = TRUE,
                    replacement = "IL\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(IL[0-9]+)(r)", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "Cd([0-9])", perl = TRUE,
                    replacement = "CD\\1", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(Ig)([a-z])", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(H[0-9])(k[0-9]+)", perl = TRUE,
                    replacement = "\\1\\U\\2", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "(B)(i+) ", perl = TRUE,
                    replacement = "\\1\\U\\2 ", ignore.case = FALSE)
  
  return(StringVec)
}

# Convenience wrapper for creating pathway enrichment plots
doPathEnrichPlot <- function(genesUp, genesDn, TERM2GENE, topN = 6) {
  # genesUp <- posGenes
  # genesDn <- negGenes
  # topN = 6
  
  EGMT <- enricher(gene = genesUp, TERM2GENE = TERM2GENE)
  eresUp <- as.data.frame(EGMT)
  if (! length(eresUp$ID)) {
    EGMT <- enricher(gene = genesUp, TERM2GENE = TERM2GENE, pvalueCutoff = .3)
    eresUp <- as.data.frame(EGMT)
  }
  EGMT <- enricher(gene = genesDn, TERM2GENE = TERM2GENE)
  eresDn <- as.data.frame(EGMT)
  if (! length(eresDn$ID)) {
    EGMT <- enricher(gene = genesDn, TERM2GENE = TERM2GENE, pvalueCutoff = .3)
    eresDn <- as.data.frame(EGMT)
  }
  eresDnIDs <- eresDn$ID
  eresDn <- eresDn[! eresDn$ID %in% eresUp$ID,]
  eresUp <- eresUp[! eresUp$ID %in% eresDnIDs,]
  eresUp <- eresUp %>% top_n(n = topN, wt = -p.adjust)
  eresUp$Group <- "Over-expressed"
  eresUp <- eresUp[order(eresUp$p.adjust),]
  eresDn$Group <- "Under-expressed"
  eresDn <- eresDn %>% top_n(n = topN, wt = -p.adjust)
  eresDn <- eresDn[order(eresDn$p.adjust),]
  pltDF <- rbind(eresUp[,c(1, 6, 10)], eresDn[,c(1, 6, 10)])
  pltDF$p.adjust <- -log10(pltDF$p.adjust)
  pltDF$EnrichmentScore <- pltDF$p.adjust
  pltDF$ID <- correlationAnalyzeR::fixStrings(pltDF$ID)
  pltDF$ID[nchar(pltDF$ID) > 40] <- paste0(substr(pltDF$ID[nchar(pltDF$ID) > 40], 1, 37), "...")
  pltDF <- pltDF[! duplicated(pltDF$ID),]
  pltDF$ID <- factor(pltDF$ID, levels = rev(pltDF$ID))
  return(ggplot(data = pltDF, 
                mapping = aes_string(x = "ID", 
                                     y = "EnrichmentScore", 
                                     fill = "Group")) +
           geom_bar(stat = "identity") + 
           ylab("-log10(pAdj)") +
           ggpubr::rotate() +
           theme_pubr(border = T, base_size = 22) +
           guides(colour = guide_legend(override.aes = list(size=3))) + 
           theme(legend.position="right") + rremove("ylab") +
           scale_fill_manual(values=c("#CC0909", "#C4BFBF")) +
           theme(axis.text.y = element_text(size = 18),
                 panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
           scale_y_continuous(expand = c(0,0), limits = c(0,(max(pltDF$EnrichmentScore)*1.1))) +
           rremove("legend"))
}

# Convenience wrapper for GSEA
myGSEA <- function(ranks,
                   TERM2GENE,
                   padjustedCutoff = .05,
                   returnDataOnly = TRUE,
                   nperm = 2000,
                   topPlots = FALSE,
                   outDir,
                   Condition = "GSEA Results",
                   plotFile = "GSEA_results") {
  
  
  # # Bug testing
  # padjustedCutoff = .05
  # topPlots = FALSE
  # nperm = 2000
  # corrDF <- correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = c("ATM"),
  #                                                   returnDataOnly = TRUE,
  #                                                   runGSEA = FALSE,
  #                                                   Sample_Type = "normal")
  # ranks <- corrDF$correlations[,1]
  # names(ranks) <- rownames(corrDF$correlations)
  # TERM2GENE <- correlationAnalyzeR::getTERM2GENE(GSEA_Type = "simple",
  #                                                Species = "hsapiens")
  
  resList <- list()
  ranks <- ranks[which(! duplicated(names(ranks)))]
  ranks <- ranks[which(! is.na(ranks))]
  ranks <- ranks[order(ranks, decreasing = TRUE)]
  
  EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                nperm = nperm, pvalueCutoff = padjustedCutoff)
  
  resGSEA <- as.data.frame(EGMT)
  
  resList[["EGMT"]] <- EGMT
  
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                  nperm = nperm, pvalueCutoff = padjustedCutoff + .15)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  if (length(resGSEA$ID) < 10){
    warning(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                   padjustedCutoff, ". Rerunning with higher pValue."))
    EGMT <- GSEA2(TERM2GENE = TERM2GENE, ranks = ranks, nproc = 1,
                  nperm = nperm, pvalueCutoff = padjustedCutoff + .45)
    resGSEA <- as.data.frame(EGMT)
    resList[["EGMT"]] <- EGMT
  }
  if (length(resGSEA$ID) < 10){
    stop(paste0("GSEA Failed -- No significant pathways at designated pValue: ",
                padjustedCutoff, ". Please check your data. If you believe this ",
                "behavior is a bug, please contact the package maintainer."))
  }
  if (topPlots) {
    resGSEA <- resGSEA[order(resGSEA$NES, decreasing = TRUE),]
    topUP <- resGSEA$ID[1:10]
    resGSEA <- resGSEA[order(resGSEA$NES, decreasing = FALSE),]
    topDOWN <- resGSEA$ID[1:10]
    resGSEA <- resGSEA[order(resGSEA$pvalue),]
    plUP <- list()
    plDOWN <- list()
    for ( i in 1:6 ) {
      pathway <- topUP[i]
      if (nchar(pathway) > 35) {
        pathTitle <- paste0(substr(pathway, 1, 30), "...")
      } else {
        pathTitle <- pathway
      }
      gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
      gp <- gp + ggplot2::labs(title = pathTitle,
                               subtitle = paste0("Enrichment score: ",
                                                 round(resGSEA$NES[which(
                                                   resGSEA$ID == pathway
                                                 )], 3)))
      gg <- ggplot2::theme_classic()
      gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                                 plot.subtitle = gg[["plot.subtitle"]],
                                 plot.margin = gg[["plot.margin"]])
      if ( i == 1) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
      } else if (i == 2) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
      } else if (i == 3) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
      } else if (i == 4) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
      } else if (i == 5) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
      } else if (i == 6) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
      }
      plUP[[i]] <- gp
      
      pathway <- topDOWN[i]
      if (nchar(pathway) > 35) {
        pathTitle <- paste0(substr(pathway, 1, 30), "...")
      } else {
        pathTitle <- pathway
      }
      gp <- clusterProfiler::gseaplot(EGMT, pathway, title = NULL)
      gp <- gp + ggplot2::labs(title = pathTitle,
                               subtitle = paste0("Enrichment score: ",
                                                 round(resGSEA$NES[which(
                                                   resGSEA$ID == pathway
                                                 )], 3)))
      gg <- ggplot2::theme_classic()
      gp <-  gp + ggplot2::theme(plot.title = gg[["plot.title"]],
                                 plot.subtitle = gg[["plot.subtitle"]],
                                 plot.margin = gg[["plot.margin"]])
      if (i == 1) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 45))
      } else if (i == 2) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 20, 20, 20))
      } else if (i == 3) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(45, 45, 20, 20))
      } else if (i == 4) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 45))
      } else if (i == 5) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 45, 20))
      } else if (i == 6) {
        gp <- gp + ggplot2::theme(plot.margin = ggplot2::margin(20, 45, 45, 20))
      }
      plDOWN[[i]] <- gp
      
    }
    
    gaUP <- ggpubr::ggarrange(plotlist = plUP, nrow = 2, ncol = 3)
    gaUP <- ggpubr::annotate_figure(gaUP,
                                    top = ggpubr::text_grob(paste0(
                                      "Top Over-Expressed Pathways in ", Condition
                                    ),
                                    size = 35)
    )
    resList[["GSEA_up"]] <- gaUP
    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gaUP,
                      filename = file.path(outDir, paste0(plotFile, "_topPathwaysUP.png")),
                      height = 14, width = 20)
    }
    
    
    gaDOWN <- ggpubr::ggarrange(plotlist = plDOWN, nrow = 2, ncol = 3)
    gaDOWN <- ggpubr::annotate_figure(gaDOWN,
                                      top = ggpubr::text_grob(paste0(
                                        "Top Under-Expressed Pathways in ", Condition
                                      ),
                                      size = 35)
    )
    
    if (! returnDataOnly) {
      ggplot2::ggsave(plot = gaDOWN,
                      filename = file.path(outDir, paste0(plotFile, "_topPathwaysDOWN.png")),
                      height = 14, width = 20)
    }
    resList[["GSEA_down"]] <- gaDOWN
    
  }
  
  if (! returnDataOnly) {
    data.table::fwrite(x = resGSEA, file = file.path(outDir,
                                                     paste0(plotFile,
                                                            "_GSEA.csv")))
  }
  
  resList[["eres"]] <- resGSEA
  
  cat("\nReturning ... ", names(resList), "\n")
  
  return(resList)
}



# Modified GSEA from clusterProfiler. Reduces computational time.
GSEA2 <- function(TERM2GENE, ranks,
                  nperm = 2000, nproc = "auto",
                  pvalueCutoff = .05) {
  if (nproc == "auto") {
    nproc = parallel::detectCores()
  }
  TERMList <- TERM2GENE %>% split(x = .$gene_symbol, f = .$gs_name)
  EGMT <- fgsea::fgsea(pathways = TERMList, nproc = nproc,
                       maxSize = 500,
                       minSize = 15,
                       stats = ranks, nperm = nperm)
  res <- data.frame(
    ID = as.character(EGMT$pathway),
    Description = as.character(EGMT$pathway),
    setSize = EGMT$size,
    enrichmentScore = EGMT$ES,
    NES = EGMT$NES,
    pvalue = EGMT$pval,
    p.adjust = EGMT$padj,
    core_enrichment = vapply(EGMT$leadingEdge, FUN.VALUE = "char",
                             paste0, collapse='/'),
    stringsAsFactors = FALSE
  )
  res <- res[!is.na(res$pvalue),]
  res <- res[ res$pvalue <= pvalueCutoff, ]
  res <- res[ res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  params <- list(pvalueCutoff = pvalueCutoff,
                 nPerm = nperm,
                 pAdjustMethod = "BH",
                 exponent = 1,
                 minGSSize = 15,
                 maxGSSize = 500
  )
  row.names(res) <- res$ID
  EGMT <- new("gseaResult",
              result     = res,
              geneSets   = TERMList,
              geneList   = ranks,
              params     = params,
              readable   = FALSE
  )
  EGMT@organism <- "UNKNOWN"
  EGMT@setType <- "UNKNOWN"
  EGMT@keytype <- "UNKNOWN"
  return(EGMT)
}
