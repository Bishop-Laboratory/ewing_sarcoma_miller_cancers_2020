######################################################################################
################################## Helper Functions ##################################
######################################################################################

# Basic function to convert mouse to human gene names (from https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/)
convertMouse2HumanGenes <- function(srtMouse, srtHuman) {
  
  # srtMouse <- GSE130146
  # srtHuman <- normalCellsInt
  
  smMouse <- srtMouse@assays$RNA@counts
  genesMouse <- rownames(smMouse)
  smHuman <- srtHuman@assays$RNA@counts
  genesHuman <- rownames(smHuman)
  if (! file.exists("Data/mouse2HumanGenes.rda")) {
    message("Gene conversion object not found -- querying BiomaRt...")
    require("biomaRt")
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    mouse2HumanGenes <- getLDS(attributes = c("mgi_symbol"), 
                      mart = mouse, attributesL = c("hgnc_symbol"),
                      martL = human, uniqueRows=T)
    save(mouse2HumanGenes, file = "Data/mouse2HumanGenes.rda")
  } else {
    load("Data/mouse2HumanGenes.rda")
  }
  mouse2HumanGenes <- mouse2HumanGenes[which(! duplicated(mouse2HumanGenes$MGI.symbol)),]
  mouse2HumanGenes <- mouse2HumanGenes[which(! duplicated(mouse2HumanGenes$HGNC.symbol)),]
  mouse2HumanGenes <- mouse2HumanGenes[mouse2HumanGenes$MGI.symbol %in% genesMouse,]
  mouse2HumanGenes <- mouse2HumanGenes[mouse2HumanGenes$HGNC.symbol %in% genesHuman,]
  genesMouse <- genesMouse[genesMouse %in% mouse2HumanGenes$MGI.symbol]
  smMouse <- smMouse[genesMouse,]
  mouse2HumanGenes2 <- mouse2HumanGenes[order(match(mouse2HumanGenes$MGI.symbol, genesMouse)),]
  all(mouse2HumanGenes2$MGI.symbol == genesMouse)
  rownames(smMouse) <- mouse2HumanGenes2$HGNC.symbol
  genesHuman <- genesHuman[genesHuman %in% mouse2HumanGenes$HGNC.symbol]
  smHuman <- smHuman[genesHuman,]
  
  srtMouse <- CreateSeuratObject(smMouse, meta.data = srtMouse@meta.data)  
  srtHuman <- CreateSeuratObject(smHuman, meta.data = srtHuman@meta.data)  
  resList <- list(
    "srtMouse" = srtMouse,
    "srtHuman" = srtHuman
  )
  return(resList)
}

# Pretty timestamp
timestamp2 <- function() {
  # Partially from https://stackoverflow.com/questions/1962278/dealing-with-timestamps-in-r
  now <- Sys.time()
  timeList <- unclass(as.POSIXlt(now))
  secNow <- ifelse(round(timeList$sec) < 10, paste0(0, round(timeList$sec)), round(timeList$sec))
  minNow <- ifelse(round(timeList$min) < 10, paste0(0, round(timeList$min)), round(timeList$min))
  paste0("[", timeList$hour,":",minNow,":",secNow,
         " ", (timeList$mon+1), "/", timeList$mday, "/", (timeList$year + 1900), "]", sep = "")
}

# Convenient wrapper for GSEA and cor testing
doTrendyAnalysis <- function(countsNorm, timeVec) {
  
  # timeVec <- timeNow
  corRes <- c()
  for (i in 1:length(countsNorm[,1])) {
    countsNow <- countsNorm[i,]
    corRes <- c(corRes, cor(x = countsNow, y = timeVec))
  }
  names(corRes) <- rownames(countsNorm)
  r2 <- corRes^2*sign(corRes)
  hist(r2, breaks = 100)
  r2 <- r2[order(r2, decreasing = FALSE)]
  # pVec <- dt(abs(corRes)/sqrt((1-corRes^2)/(n-2)), df = 2)
  pltDF <- data.frame(
    "geneName" = names(r2),
    "value" = r2,
    "rank" = seq(1, length(r2)),
    stringsAsFactors = F
  )
  pltDF <- pltDF[order(pltDF$rank, decreasing = T),]
  return(pltDF)
  # GSEARes <- myGSEA(r2, TERM2GENE, padjustedCutoff = pVal)
  # GSEARes[["rankDF"]] <- pltDF
  # return(GSEARes)
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
featurePlot2 <- function(pltData, countsNorm, features, featureName = NULL,
                         plotName,  ptSize = .4, normalize = FALSE,
                         height = 7.5, width = 8, legend = TRUE,
                         x, y, logNorm = FALSE, doTopCut = FALSE,
                         pal = colorRampPalette(brewer.pal(9, "OrRd"))(255)) {
  # pltData <- pltData
  # countsNorm <- vsd
  # x <- "UMAP_1"
  # y <- "UMAP_1"
  # normalize <- FALSE
  # doTopCut = F
  # logNorm = F
  # ptSize = .4
  # pal = colorRampPalette(brewer.pal(9, "OrRd"))(255)
  # features <- bachMarkers
  # signature <- upLoops
  # pal = colorRampPalette(brewer.pal(9, "OrRd"))(255)
  
  if (! "samples" %in% colnames(pltData)) {
    if (! rownames(pltData) == colnames(countsNorm)) {
      stop("Check to make sure the 'samples' column exists in pltData")
    } else {
      pltData$samples <- rownames(pltData)
    }
  }
  countsNorm <- countsNorm[,colnames(countsNorm) %in% pltData$samples]
  if (dim(countsNorm)[2] == 0) {
    stop("Check to make sure the 'samples' column in pltData is equal to colnames of countData")
  }
  pltData <- pltData[pltData$samples %in% colnames(countsNorm),]
  pltData <- pltData[order(match(pltData$samples, colnames(countsNorm))),]
  if (normalize) {
    Sizes <- MedianNorm(countsNorm)
    countsNorm <- GetNormalizedMat(countsNorm, Sizes)
  }
  all(pltData$samples == colnames(countsNorm))
  features <- features[features %in% rownames(countsNorm)]
  pltData$plotFeat <- colMeans(countsNorm[features,, drop = F])
  colorName <- ifelse(! is.null(featureName), featureName, ifelse(
    length(features) > 1, "Signature\nexpression", paste0(features[1], "\nexpression")
  ))
  pltOrd <- pltData$plotFeat
  pltOrd[is.na(pltOrd)] <- 0
  pltData <- pltData[order(pltOrd, decreasing = FALSE),]
  g1 <- ggplot(pltData, mapping = aes_string(x = x, y = y,
                                                color = pltData$plotFeat)) +
    geom_point(size = ptSize) + theme_pubr(border = T, base_size = 22) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    # scale_color_viridis(option = "B") +
    theme(legend.position="right") +
    rremove("legend") +
    labs(color = colorName) +
    scale_color_gradientn(colors = pal, na.value = "#DBDBDB")
  
  ggsave(g1, filename = paste0("Figures_v2/", plotName, ".png"), height = height, 
         width = width)
  
  if (legend) {
    g1 <- ggplot(pltData, mapping = aes_string(x = x, y = y,
                                                  color = pltData$plotFeat)) +
      geom_point(size = ptSize) + theme_pubr(border = T, base_size = 22) +
      theme(#border = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
      # scale_color_viridis(option = "B") +
      theme(legend.position="right") +
      # rremove("legend") +
      labs(color = colorName) +
      scale_color_gradientn(colors = pal, na.value = "#DBDBDB")
    ggsave(g1, filename = paste0("Figures_v2/", plotName, "_legend.png"), height = height, 
           width = round((width * 1.33), digits = 2))
  }
  return(g1)
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
                     ptSize = .4, colorMap = NULL, xlab = NULL, ylab = NULL,
                     legend = TRUE, silent = TRUE) {
  g1 <- ggplot(data = data, mapping = mapping) + geom_point(size = ptSize) +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
  if (! is.null(colorMap)) {
    g1 <- g1 +
      scale_color_manual(values = colorMap)
  }
  if (! is.null(xlab)) {
    g1 <- g1 + xlab(xlab)
  }
  if (! is.null(ylab)) {
    g1 <- g1 + ylab(ylab)
  }
  ggsave(g1, filename = paste0("Figures_v2/", plotName, ".png"), height = height, 
         width = width)
  if (legend) {
    g1 <- ggplot(data = data, mapping = mapping) + geom_point(size = ptSize) + theme_pubr(border = T, base_size = 22) +
      guides(colour = guide_legend(override.aes = list(size=3))) + 
      theme(legend.position="right") +
      theme(#border = element_line(colour = "black", size = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
    if (! is.null(colorMap)) {
      g1 <- g1 +
        scale_color_manual(values = colorMap)
    }
    if (! is.null(xlab)) {
      g1 <- g1 + xlab(xlab)
    }
    if (! is.null(ylab)) {
      g1 <- g1 + ylab(ylab)
    }
    ggsave(g1, filename = paste0("Figures_v2/", plotName, "_legend.png"), height = height, 
           width = round((width * 1.33), digits = 2))
  }
  if (! silent) {
    return(g1)
  } 
}

# Convenience function for PHATE correlation plots
doPhateCorrPlot <- function(nameNow, resSamples,
                            intGenes = c("FANCI", "BRCA1", "CD44")) {
  # resSamples <- resListEWS
  # nameNow <- "PHATE_1"
  
  phateRes <- resSamples[[nameNow]]
  pltDF <- phateRes
  pltDF$gene <- ""
  m <- length(pltDF$rank)-4
  n <- max(pltDF$value)
  o <- min(pltDF$value)
  pltDF$gene[(pltDF$geneName %in% intGenes |
                pltDF$rank < 4 | pltDF$rank > m)] <-as.character(pltDF$geneName)[(pltDF$geneName %in% intGenes|
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

doGSEABarPlot <- function(eres, topN = 5, charLimit = 50, forcePval = FALSE,
                          colors = c("firebrick", "skyblue")) {
  if (! forcePval) {
    eres <- eres[eres$p.adjust < .05,]
    if (! length(eres$ID)) {
      warning("No significant pathways returned")
      return(NULL)
    }
  }
  
  eresUp <- eres %>% top_n(n = topN, wt = NES)
  eresUp <- eresUp[eresUp$NES > 0,]
  eresUp$Group <- "Over-expressed"
  eresDn <- eres %>% top_n(n = topN, wt = -NES)
  eresDn <- eresDn[eresDn$NES < 0,]
  eresDn$Group <- "Under-expressed"
  eresPlt <- rbind(eresUp, eresDn)
  eresPlt <- eresPlt[order(eresPlt$NES, decreasing = T),]
  eresPlt$ID <- fixStrings(eresPlt$ID)
  eresPlt$ID[nchar(eresPlt$ID) > charLimit] <- paste0(substr(eresPlt$ID[nchar(eresPlt$ID) > charLimit], 1, (charLimit-3)), "...")
  eresPlt <- eresPlt[! duplicated(eresPlt$ID),]
  eresPlt$ID <- factor(eresPlt$ID, levels = rev(eresPlt$ID))
  g1 <- ggplot(data = eresPlt, 
               mapping = aes_string(x = "ID", 
                                    y = "NES", 
                                    fill = "Group")) +
    geom_bar(stat = "identity") + 
    ylab("Normalized Enrichment") +
    ggpubr::rotate() +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") + rremove("ylab") +
    scale_fill_manual(values=colors) +
    theme(axis.text.y = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_y_continuous(expand = c(0,0), limits = c((min(eresPlt$NES)*1.1),(max(eresPlt$NES)*1.1))) +
    rremove("legend")
  return(g1)
}

doAltGSEA <- function(ranks, TERM2GENEList, nCharLimit = 50, nGS = 5, pCut = .05) {
  res <- lapply(names(TERM2GENEList), FUN = function(TERM2GENENowName) {
    TERM2GENENow <- TERM2GENEList[[TERM2GENENowName]]
    GSEARes <- myGSEA(ranks = ranks, TERM2GENE = TERM2GENENow, padjustedCutoff = pCut)
    list("gseaPlot" = (doGSEABarPlot(GSEARes$eres, charLimit = nCharLimit, topN = nGS, colors = c("#525252", "#BDBDBD")) +
                         labs(title = TERM2GENENowName)),
         "eres" = GSEARes$eres)
  })
  names(res) <- names(TERM2GENEList)
  return(res)  
}

# Taken from https://github.com/cran/VennDiagram/blob/master/R/hypergeometric.test.R
#This function performs the hypergeometric test on the two categories. Taken from package BoutrosLab.statistics.general
calculate.overlap.and.pvalue = function(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) {
  
  # calculate actual overlap
  actual.overlap <- length(intersect(list1, list2));
  
  # calculate expected overlap
  # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
  expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;
  
  adjust.value <- 0;
  
  # adjust actual.overlap to reflect P[X >= x]
  if (adjust & !lower.tail) {
    adjust.value <- 1;
    warning('Calculating P[X >= x]');
  }
  
  # calculate significance of the overlap
  overlap.pvalue <- phyper(
    q = actual.overlap - adjust.value,
    m = length(list1),
    n = total.size - length(list1),
    k = length(list2),
    lower.tail = lower.tail
  );
  
  # return values
  return( c(actual.overlap, expected.overlap, overlap.pvalue) );
  
}



# Convenience wrapper for msigdbr
getTERM2GENE <- function(GSEA_Type = c("simple"),
                         Species = c("hsapiens", "mmusculus"),
                         sampler = FALSE) {
  # require(dplyr)
  # require(tidyr)
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
    dplyr::filter(.data$gs_cat %in% categories) %>%
    dplyr::select(.data$gs_name, .data$gene_symbol)
  
  if (sampler) {
    print("Using sampler!")
    set.seed(1)
    TERM2GENE <- TERM2GENE[sample(nrow(TERM2GENE), size = 100000),]
  }
  return(TERM2GENE)
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
  # StringVec <- "IGLESIAS_E2F_TARGETS_UP"
  
  # These SILIGAN genesets are inverted... 
  StringVec <- gsub(StringVec, pattern = "SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN",
                    replacement = "SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_DN_(UP)")
  StringVec <- gsub(StringVec, pattern = "SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_UP",
                    replacement = "SILIGAN_TARGETS_OF_EWS_FLI1_FUSION_UP_(DN)")
  # So are the IGLESIAS... 
  StringVec <- gsub(StringVec, pattern = "IGLESIAS_E2F_TARGETS_UP",
                    replacement = "IGLESIAS_E2F_TARGETS_UP_(DN)")
  StringVec <- gsub(StringVec, pattern = "IGLESIAS_E2F_TARGETS_DN",
                    replacement = "IGLESIAS_E2F_TARGETS_UP_(UP)")
  
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
  StringVec <- gsub(StringVec, pattern = " Dn \\(", replacement = " Down (", ignore.case = FALSE)
  StringVec <- gsub(StringVec, pattern = "\\(Dn\\)", replacement = "(Down)", ignore.case = FALSE)
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
  StringVec <- gsub(StringVec, pattern = "IgLesias", perl = TRUE,
                    replacement = "Iglesias", ignore.case = FALSE)
  return(StringVec)
}

# Convenience wrapper for creating pathway enrichment plots
doPathEnrichPlot <- function(genesUp, genesDn, TERM2GENE, returnData = FALSE, ncharLimit = 50,
                             topN = 6, colors = c("#CC0909", "#C4BFBF")) {
  # genesUp <- posGenes
  # genesDn <- negGenes
  # topN = 6
  
  EGMT <- enricher(gene = genesUp, TERM2GENE = TERM2GENE)
  eresUpRaw <- as.data.frame(EGMT)
  if (! length(eresUpRaw$ID)) {
    EGMT <- enricher(gene = genesUp, TERM2GENE = TERM2GENE, pvalueCutoff = .3)
    eresUpRaw <- as.data.frame(EGMT)
  }
  EGMT <- enricher(gene = genesDn, TERM2GENE = TERM2GENE)
  eresDnRaw <- as.data.frame(EGMT)
  if (! length(eresDnRaw$ID)) {
    EGMT <- enricher(gene = genesDn, TERM2GENE = TERM2GENE, pvalueCutoff = .3)
    eresDnRaw <- as.data.frame(EGMT)
  }
  eresDnIDs <- eresDnRaw$ID
  eresDn <- eresDnRaw[! eresDnRaw$ID %in% eresUpRaw$ID,]
  eresUp <- eresUpRaw[! eresUpRaw$ID %in% eresDnIDs,]
  eresUp <- eresUp %>% top_n(n = topN, wt = -p.adjust)
  eresUp$Group <- "Over-expressed"
  eresUp <- eresUp[order(eresUp$p.adjust),]
  eresDn$Group <- "Under-expressed"
  eresDn <- eresDn %>% top_n(n = topN, wt = -p.adjust)
  eresDn <- eresDn[order(eresDn$p.adjust),]
  pltDF <- rbind(eresUp[,c(1, 6, 10)], eresDn[,c(1, 6, 10)])
  pltDF$p.adjust <- -log10(pltDF$p.adjust)
  pltDF$EnrichmentScore <- pltDF$p.adjust
  pltDF$ID <- fixStrings(pltDF$ID)
  pltDF$ID[nchar(pltDF$ID) > ncharLimit] <- paste0(substr(pltDF$ID[nchar(pltDF$ID) > ncharLimit], 1, (ncharLimit - 3)), "...")
  pltDF <- pltDF[! duplicated(pltDF$ID),]
  pltDF$ID <- factor(pltDF$ID, levels = rev(pltDF$ID))
  g1 <- ggplot(data = pltDF, 
               mapping = aes_string(x = "ID", 
                                    y = "EnrichmentScore", 
                                    fill = "Group")) +
    geom_bar(stat = "identity") + 
    ylab("-log10(pAdj)") +
    ggpubr::rotate() +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") + rremove("ylab") +
    scale_fill_manual(values=colors) +
    theme(axis.text.y = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,(max(pltDF$EnrichmentScore)*1.1))) +
    rremove("legend")
  if (! returnData) {
    return(g1)
  } else {
    return(list("plot" = g1,
                "eresUp" = eresUpRaw,
                "eresDn" = eresDnRaw))
  }
  
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
  set.seed(42)
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


# Convenience wrapper for making supplemental cluster marker plots
doClusterMarkerPlot <- function(srt, srtMarkers, outName, ncharLimit = 50,
                                intClust = NULL, topN = 6,
                                TERM2GENE = NULL) {
  # # Bug testing
  # srt <- BDInt
  # srtMarkers <- BDIntMarkers
  # outName = "Figures_v2/Fig1Cx0iii"
  
  if (is.null(TERM2GENE)) {
    TERM2GENE <- getTERM2GENE(GSEA_Type = "simple")
  }
  
  dir.create(paste0("Figures_v2/", outName), showWarnings = F)
  srtMarkers <- srtMarkers[srtMarkers$p_val_adj < .05,]
  srtMarkers$p_val_adj[srtMarkers$p_val_adj == 0] <- .Machine$double.xmin
  clusterList <- unique(as.numeric(srt$seurat_clusters))
  clusterList <- c(0, clusterList[order(clusterList)])
  names(clusterList) <- paste0("Cluster: ", clusterList)
  pltDat <- srt@meta.data
  pltDat$UMAP_1 <- srt@reductions$umap@cell.embeddings[,c(1)]
  pltDat$UMAP_2 <- srt@reductions$umap@cell.embeddings[,c(2)]
  for (i in 1:length(clusterList)) {
    if (! is.null(intClust)) {
      if (! i %in% (intClust+1)) {
        next
      }
    }
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
              fixed.asp = T, rot.per = 0, mar = c(0,0,0,0),
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
    g4 <- doPathEnrichPlot(posGenes, negGenes, topN = topN, ncharLimit = ncharLimit, 
                           TERM2GENE = TERM2GENE)
    # Arrange
    ggNow <- ggdraw() +
      draw_plot(g1, x = 0, y = .4, width = .4, height = .6) +
      draw_plot(g2, x = 0, y = 0, width = .2, height = .35) +
      draw_plot(g3, x = 0.25, y = 0, width = .2, height = .35) +
      draw_plot(g4, x = 0.47, y = 0, width = .53, height = 1) 
    ggsave(ggNow, file = paste0("Figures_v2/", outName, "/Cluster_", clusterNow, ".png"),
           height = 9, width = 16)
  }
}


# Tissue dict
tissueDict <- list("brain" = list("yes" = c("cortex", "brain", "lobe", "hippoc", "^pfc$", "mge",
                                            "^vc$", "^cbc$", "gyrus", "stroke", "sciencell",
                                            "alzheim", "frontal", "dentate", "white matter", "brian",
                                            "cranial",
                                            "grey matter", "gray matter", "striatum", "pericyte",
                                            "nerv", "gangli", "bipol", "medull", "putamen",
                                            "hippocamp", "neur", "glia", "amygdala", "oligodendro",
                                            "spine", "spinal", "astrocyt", "cereb"),
                                  "no" = c("liver", "kidney", "microgli", "aneurysm", "vessel",
                                           "precursor", 'progenitor', "stem cell", "Neuroectodermal",
                                           "iPSC", "NPC", "NSC")),
                   "thyroid" = c("thyroid"),
                   "respiratory" = c("lung", "airway", "nasal", "hsaec",
                                     "trach", "pleura", "alveol", "bronch"),
                   "skin" = c("skin", "keratin", "dermis", "^dk$", "epidermis", "melano", "psoriasis"),
                   "pancreas" = list("yes" = c("pancreas", "pancrea", 
                                               "islet", "alpha", "beta", "delta", "epsilon"),
                                     "no" = c("falpha", "fbeta")),
                   "kidney" = list("yes" = c("kidney", "nephr", "glomerul", 
                                             "renal", "clear cell", "^ptec$"),
                                   "no" = c("ASDLJSND")),
                   "fetal" = c("placent", "fetal", "huvec", "decidua", "germ",
                               "fetus", "embry", "umbil", "cord blood"),
                   "cartilage" = c("cartilag", "chondr", "joint"),
                   "mammary" = list("yes" = c("mammary", "breast", "imec", "hmec",
                                              "reduction mammoplasty no known cancer"),
                                    "no" = c("ASDLKNASDL")),
                   "stomach"= c("stomach", "gastric"),
                   "esophagus" = c("esophag"),
                   "intestines" = list("yes" = c("intestine", "intestinal", "colon", "duoden", "colorect",
                                                 "ileum", "ileal", "gut", "bowel", "jejun", "sigmoid",
                                                 "recto", "ileocolic"),
                                       "no" = c("colonization", "human colorectal cell line", 
                                                "ncm356d")),
                   "muscle" = list("yes" = c("muscle", "myo",  "^smc", " smc",
                                             "lateralis", "gastrocnemius",
                                             "skeletal", "brach", "satellite", "ceps"), 
                                   "no" = c("cardiac", 'heart', "endothel", "vessel")),
                   "liver" = list("yes" = c("liver", "hepat", "kupffer", "phh"),
                                  "no" = c("deliver")),
                   "adipose" = list(
                     "yes" = c("adipose", "fat", "adipo", "^wat$", "^bat$"),
                     "no" = c("mesenchymal", "milk", "adipocyterna", "hdfatprx1")
                   ),
                   "stem-like" = list("yes" = c("stem", "progen", "prog$", "blast ", 
                                                "hesc", "hues64", "blast$",
                                                "cd34", "cord blood", "wharton", "pluripotent",
                                                "ncsc", "ncc", "hnspc", "^kp$", "ECFC",
                                                "npc", "h1esc", "^h1", " pes[0-9]+", "hff",
                                                "embryo", "ectoderm", "epsc", "morula",
                                                "mesoderm", "^ips$", "hpsc", "hfl1", "HSPC",
                                                "dental pulp cells", "cpcs", "fbs", "primitive",
                                                "[a-zA-Z]genic", "hematopoietic", "^esc$",
                                                "human es", "human ips", "endoderm", "sscs",
                                                "human icm", "human te", "derma", "^esc$",
                                                "NAMEC", "mesenchym", "hffs", "hmsc",
                                                "h9", "controlc7", "stroma", "msc", "zygote",
                                                "fibroblas", "IMR90", "hff1", "nhdf",
                                                "npsc", "ipsc", "poetic", "ips cell", "satellite"),
                                      "no" = c("dermal")),
                   "cardiac" = list("yes" = c("cardiac", "heart", 
                                              "atria", "atrium",
                                              "coron", "aort", "ventric"),
                                    "no" = c("SADONASDIOFN")),
                   "endothelial" = list("yes" = c("endoth", "huvec", "hdmec", "ECFC",
                                                  "vascul", "vessel","f[0-9]ecs",
                                                  "hpmec", "ecctr", "ectnf", "ecil"),
                                        "no" = c("ASDKFBASD")),
                   "spleen" = c("spleen", "splen"),
                   "bladder" = list("yes" = c("bladder", "urin", "urothe"),
                                    "no" = c("during")),
                   "retina" = c("retina", "macular", "retin", "photo"),
                   "thymus" = c("thymus", "thymic"),
                   "male reproductive" = list("yes" = c("testis", "testes", 
                                                        "leydig", "peritubular", "sertoli",
                                                        "cauda", "corpus", "caput",
                                                        "testicle", "sperm",
                                                        "epidid", "gonad"),
                                              "no" = c("prostate", "prostatic")),
                   "prostate" = list("yes" = c("prostate", "prostatic"),
                                     "no" = c("ASODINASD")),
                   "female reproductive" = list("yes" = c("ovar",  "amniotic", 
                                                          "uter", "placent", "fallopian",
                                                          "deciduo", "cervix", "amnion", "chorionic villus",
                                                          "endometri", "oviductal",
                                                          "cervic", "vagi", "granulosa"),
                                                "no" = c("stroma")),
                   "immune" = list("yes" = c("immune", "macroph", "leuk",
                                             "killer", "lymph",
                                             "cd[0-9]+",
                                             "NKT", "blood", "microgli", "gm12878", "tcell", "bcell",
                                             "b cell", "t cell", "cd4", "cd8", "nk cell",
                                             "monocyt", "dendrit", "granulocyt",
                                             "lympho", "mononucle", "pbmc", "neutro", "treg"),
                                   "no" = c("ctc", "stroma", "osteo", "vein",
                                            "vessel", "cord", "oma", "cd34",
                                            "[a-zA-Z]b cell", "[a-zA-Z]t cell")),
                   "bone" = list("yes" = c("femur", "osteo", "hfob",
                                           "mandible", "bone", "joint"),
                                 "no" = c("MSC")),
                   "ewing sarcoma" = list("yes" = c("ewing", "a673", "$tc[0-9]+",
                                                    " tc[0-9]+", "ews", "ews502",
                                                    "chla[0-9]+", "TC32",
                                                    "tc71", "sknmc", 
                                                    "rd-es", "sk-es", 
                                                    "skes", "sk-nm-c",
                                                    "cadoes", "^rdes$",
                                                    "^rdes ",
                                                    "^rdes_", " es2 ",
                                                    "^es2 ", " es2$",
                                                    "ew[0-9]+", "es7"),
                                          "no" = c("pES", "RUES", "HUES", 
                                                   "PES1", "hes", "CSES7", 
                                                   "A375", "786", "TTC1240", "TTC642",
                                                   "hela", "HepG2", "k562", "RWPE1", "aES7")),
                   "neural crest cells" = list("yes" = c("ncc", "ncsc"),
                                               "no" = c("NCCIT", "NCC24")),
                   "HSCs" = list("yes" = c("HSC", "hematopoeitic stem", "cd34"),
                                 "no" = c('hepatic', "shScramble")),
                   "MSCs" = list("yes" = c("MSC", "mesenchymal", "mesenchymal"),
                                 "no" = c("UMSCC", "pMSCV")),
                   "tumor" = list("yes" = c("cancer", "carcin", "sarcom", "metasta", "tumor",
                                            "[a-zA-Z]oma ", "[a-zA-Z]oma$"),
                                  "no" = c("healthy", "normal", "stroma")),
                   "Neural progenitor/stem cells" = c("NPC", "NSC", "npsc", "pMN progenitor",
                                                      "hnspc", "Neuron Precursor", "neural precursor"),
                   "hESCs" = list("yes" = c("hesc", "^esc$", "^h9$", " h9 ", 
                                            " h9$", "^h9 ", "h1esc",
                                            "pES[0-9]+","HUES[0-9]+", "embyonic", 
                                            "human es", "^h1", "embryonic stem"),
                                  "no" = c("derived", "NPC", "MSC")),
                   "iPSCs" = list("yes" = c("iPSC", "ips cell", "^ips$","hiPS",
                                            "human ips", "hpsc", "pluripotent stem", "Induced pluripotent"),
                                  "no" = c("iPSC-derived", "NPC", "MSC", "ips derived", "ipsc derived",
                                           "ips-derived")),
                   "Other stem-cells" = c("pscs", "cpcs", 
                                          "dental pulp", "periodontal ligament stem"),
                   "Other prenatal tissues" = list("yes" = c("parthenogenic", "blastoycst", "trophoblast",
                                                             "zygote", "endoderm", "morula", "oocyte", "Umbilical",
                                                             " ICM$", " te$", "oocy", "oophorus", 
                                                             "hek293", "hek 293", " 293T$", "^293T ",
                                                             " 293$", "^293 ",
                                                             "primitive streak", "fetal", "Germinal",
                                                             "Embryo", "Fetus","mesoderm","ectoderm", "endometrial stroma"),
                                                   "no" = c("cord blood")),
                   "fibroblasts" = c("hff1", "HFF-1", "nhdf", "^derma$",
                                     "fibroblas", "IMR90", "HFF", "fbs"),
                   "Non-Ewing tumor" = list("yes" = c("cancer", "[a-z]oma", "hela", "metast",
                                                      "mdamd231", "^rt[0-9]+", "g401", "Soft Tissue, Mesenchymal",
                                                      "tumor", "hec1b", "^omental tissue$", "HNSCC", "HCC",
                                                      "hela", "k562", "reh", "jurkat", "leukemi", "293", "bewo",
                                                      "kras", "mcf", "lncap", "bjab", "gbm", " aml",
                                                      "rko", "ramos", "mel888", "aml ",
                                                      "vcap", "saos2", "vapc", "nalm6", "set2", "tov21",
                                                      "cancer", "carcin", "sarcom", "metasta", "tumor",
                                                      "[a-zA-Z]oma ", "[a-zA-Z]oma$", "NCCIT",
                                                      "panc1", "mcf7", "pc3"),
                                            "no" = c("healthy", "normal","ewing", "a673", 
                                                     "tc32", "ews", "ews502", "es2",
                                                     "chla", "tc71", "sknmc", 
                                                     "SK-NM-C", "SK-ES",
                                                     "cadoes", "^rdes$",
                                                     "^rdes_", "stroma",
                                                     "ew[0-9]+", "es7")),
                   "SingleCell" = c("single cell", "single-cell", "smart seq", "in-drop",
                                    "cel-seq", "10X", "scRNA seq", "smartseq", "CELseq",
                                    "smart-seq", "indrop", "drop-seq", "drop seq", "single nucleus",
                                    "single-nucleus", "snRNA-Seq", "snRNASeq",
                                    "fluidigm", "scRNASeq", "scRNA-Seq", "chromium"))

jsonlite::write_json(x = tissueDict, path = "Data/bulkRNASeq/tissueDictionary.json")


