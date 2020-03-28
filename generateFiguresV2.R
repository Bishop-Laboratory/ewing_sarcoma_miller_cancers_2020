#####################################################################################
#####################################################################################
#   _____                            _         ______ _                             #
#  / ____|                          | |       |  ____(_)                            #
#  | |  __  ___ _ __   ___ _ __ __ _| |_ ___  | |__   _  __ _ _   _ _ __ ___  ___   #
#  | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \ |  __| | |/ _` | | | | '__/ _ \/ __|  #
#  | |__| |  __/ | | |  __/ | | (_| | ||  __/ | |    | | (_| | |_| | | |  __/\__ \  #
#   \_____|\___|_| |_|\___|_|  \__,_|\__\___| |_|    |_|\__, |\__,_|_|  \___||___/  #
#                                                        __/ |                      #
#                                                       |___/                       #
#####################################################################################
#####################################################################################

# R Script for fully reproducing Miller et al. 2020

# To run:
## 1. Clone github repo at: https://github.com/millerh1/Ewing-sarcoma-paper-Miller-2020.git
## 2. cd Ewing-sarcoma-paper-Miller-2020/
## 3. conda env create -f environment.yml
## 4. conda activate ewsPaperEnv
## 5. (Rscript generateFigures.R) |& tee generateFigures_logFile.txt

#####################################################################################
#####################################################################################

cat("\nReproducing analysis from Miller et al. 2020\n")
source("helpers_v2.R")

#####################################################################################
############################ Preliminary: Load libraries ############################ 
#####################################################################################

cat("\n", timestamp2(), " Loading libraries and other dependencies...\n", sep = "")
suppressMessages(library(ggpubr))
suppressMessages(library(rhdf5))
suppressMessages(library(DESeq2))
suppressMessages(library(cowplot))
suppressMessages(library(uwot))
suppressMessages(library(biomaRt))
if (! "phateR" %in% installed.packages()[,1]) {
  install.packages("phateR", repos = "https://cloud.r-project.org/")
}
suppressMessages(library(phateR))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(ggplotify))
suppressMessages(library(RcppParallel))
suppressMessages(library(dplyr))
suppressMessages(library(doParallel))
suppressMessages(library(future))
suppressMessages(library(foreach))
suppressMessages(library(scatterplot3d))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(DiffBind))
suppressMessages(library(Rmagic))
suppressMessages(library(viridis))
suppressMessages(library(plot3D))
suppressMessages(library(wordcloud))
suppressMessages(library(ggridges))
suppressMessages(library(plotly))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ChIPpeakAnno))
suppressMessages(library(EBSeq))
suppressMessages(library(ggrepel))
suppressMessages(library(VennDiagram))
suppressMessages(library(Seurat))
suppressMessages(library(jsonlite))
# Get annotation objects
TERM2GENE <- suppressWarnings(getTERM2GENE())
TERM2GENE_GO <- suppressWarnings(getTERM2GENE(GSEA_Type = "GO:BP"))
TERM2GENE_TF <- suppressWarnings(getTERM2GENE("TF targets"))
TERM2GENE_CGP <- suppressWarnings(getTERM2GENE("Perturbations"))
RIGGI2014 <- read.csv("Data/RIGGI_2014_geneSet.csv", stringsAsFactors = FALSE)
AYNAUD2020ICEWS <- read.csv("Data/Aynaud2020_IC_EwS.csv", stringsAsFactors = FALSE)
AYNAUD2020ICEMT <- read.csv("Data/AYNAUD2020_IC10_IC30.csv", stringsAsFactors = FALSE)
newEWSSets <- data.frame(
  "gs_name" = c(rep("RIGGI_2014_ACTIVATED_BY_EWS_FLI1", length(RIGGI2014$Activated)),
                rep("RIGGI_2014_REPRESSED_BY_EWS_FLI1", length(RIGGI2014$Repressed)),
                rep("AYNAUD_2020_ACTIVATED_BY_EWS_FLI1", length(AYNAUD2020ICEWS$Gene.Symbol[AYNAUD2020ICEWS$Prioritized.candidate.direct.target == 1])),
                rep("AYNAUD_2020_EWS_FLI1_HIGH_MARKERS", length(AYNAUD2020ICEWS$Gene.Symbol)),
                rep("AYNAUD_2020_EWS_FLI1_LOW_MARKERS", length(AYNAUD2020ICEMT$IC30_EMT))),
  "gene_symbol" = c( RIGGI2014$Activated, RIGGI2014$Repressed, 
                     AYNAUD2020ICEWS$Gene.Symbol[AYNAUD2020ICEWS$Prioritized.candidate.direct.target == 1],
                     AYNAUD2020ICEWS$Gene.Symbol, AYNAUD2020ICEMT$IC30_EMT)
)
TERM2GENE_EWS <- TERM2GENE_CGP[grep(TERM2GENE_CGP$gs_name, pattern = "EWING|_EWS|EWSR1"),]
TERM2GENE_EWS <- rbind(TERM2GENE_EWS, newEWSSets)
rmPaths <- c("RORIE_TARGETS_OF_EWSR1_FLI1_FUSION_DN", # Neuroblastoma is not an appropriate control for a EWS vs normal tissue study
             "RORIE_TARGETS_OF_EWSR1_FLI1_FUSION_UP",
             "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_UP", # These are in mice and are from hematopoetic progenitors
             "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_TOP20_UP",
             "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_DN",
             "TORCHIA_TARGETS_OF_EWSR1_FLI1_FUSION_TOP20_DN")
TERM2GENE_EWS <- unique(TERM2GENE_EWS[!TERM2GENE_EWS$gs_name %in% rmPaths,])
TERM2GENEList <- list("GO:BP" = TERM2GENE_GO,
                      "CGP" = TERM2GENE_CGP,
                      "Transcription Factors" = TERM2GENE_TF,
                      "Ewing sarcoma" = TERM2GENE_EWS)

# Reproducible...
set.seed(42)
dir.create("Figures_v2/", showWarnings = FALSE)
maxCores <- floor(defaultNumThreads()*.65)

#####################################################################################
########################### Preliminary: Download Datasets ########################## 
#####################################################################################

cat("\n", timestamp2(), " Step 2: Downloading missing datasets\n", sep = "")
if (! file.exists("Data/bulkRNASeq/human_matrix.h5")) {
  cat("\nARCHS4 expression data: \n")
  download.file(url = "https://mssm-seq-matrix.s3.amazonaws.com/human_matrix.h5", 
                destfile = "Data/bulkRNASeq/human_matrix.h5")
}
if (! file.exists("Data/bulkRNASeq/cellosaurus.txt")) {
  cat("\nCellosaurus ontology file: \n")
  download.file(url = "ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt", 
                destfile = "Data/bulkRNASeq/cellosaurus.txt")
}
if (! file.exists("Data/DRIPSeq/version0.6.0/genome.fa.bwt")) {
  cat("\niGenomes BWA 0.6.0 indices")
  download.file(url = "http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/NCBI/GRCh38/Homo_sapiens_NCBI_GRCh38.tar.gz", 
                destfile = "Data/DRIPSeq/iGenomeGRCh38.tar.gz")
  untar("Data/DRIPSeq/iGenomeGRCh38.tar.gz", compressed = TRUE, verbose = TRUE,
        exdir = "Data/DRIPSeq/")
  file.copy(from = "Data/DRIPSeq/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/version0.6.0", 
            to = "Data/DRIPSeq/", recursive = TRUE)
  unlink("Data/DRIPSeq/Homo_sapiens/", recursive = TRUE, force = TRUE)
  file.remove("Data/DRIPSeq/iGenomeGRCh38.tar.gz")
}
dir.create("Data/scRNASeq/CountMatrices", showWarnings = FALSE)
dir.create("Data/scRNASeq/CountMatrices/Delattre_EWS_PDX", showWarnings = FALSE)
dir.create("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells", showWarnings = FALSE)
if (! file.exists("Data/scRNASeq/CountMatrices/Delattre_EWS_PDX/GSM3730317_PDX-861_matrix.mtx.gz")) {
  download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130025/suppl/GSE130025_RAW.tar", destfile = "Data/scRNASeq/CountMatrices/Delattre_EWS_PDX.tar")
  untar("Data/scRNASeq/CountMatrices/Delattre_EWS_PDX.tar", exdir = "Data/scRNASeq/CountMatrices/Delattre_EWS_PDX")
  rmFiles <- list.files("Data/scRNASeq/CountMatrices/Delattre_EWS_PDX/", pattern = "\\.txt|\\.bed", full.names = T)
  file.remove(rmFiles)
}
if (! file.exists("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368463_CHLA10_matrix.mtx.gz")) {
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE146221&format=file", destfile = "Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells.tar")
  untar("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells.tar", exdir = "Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells")
}

#####################################################################################
#################################### Bulk RNASeq ####################################
#####################################################################################

cat("\n", timestamp2(), " Step 3: Generate Bulk RNA-Seq figures\n", sep = "")
cat("\n", timestamp2(), " Filtering unwanted samples and categorizing expression metadata...\n", 
    sep = "")

dataFile <- "Data/bulkRNASeq/human_matrix.h5"

if (! file.exists("Data/bulkRNASeq/colDataFinal.rda")) {
  # Extract metadata
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  genes <- h5read(dataFile, name = "meta/genes")
  tissue <- h5read(dataFile, name = "meta/Sample_source_name_ch1")
  description <- h5read(dataFile, name = "meta/Sample_description")
  series <- h5read(dataFile, name = "meta/Sample_series_id")
  organism <- h5read(dataFile, name = "meta/Sample_organism_ch1")
  institute <- h5read(dataFile, name = "meta/Sample_contact_institute")
  extractProtocol <- h5read(dataFile, name = "meta/Sample_extract_protocol_ch1")
  reads_aligned <- h5read(dataFile, name = "meta/reads_aligned")
  
  # Build colData
  colData <- data.frame(samples, tissue, series, organism, institute,
                        extractProtocol, reads_aligned)
  
  # Get REGEX dictionary
  tissueDictionary <- read_json("Data/bulkRNASeq/tissueDictionary.json")
  
  # Filter out single cell datasets
  patternSingle <- unlist(tissueDictionary$SingleCell)
  colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), 
                                 invert = T, perl = T,
                                 x = colData$extractProtocol, ignore.case = T)),]
  colData$tissue <- gsub(colData$tissue, pattern = "_", replacement = " ")
  
  # Categorize cancer samples
  cancerString <- unlist(tissueDictionary$tumor[[1]])
  nonCancerString <- c(unlist(tissueDictionary$tumor[[2]]), 
                       unlist(tissueDictionary$iPSCs[[1]]),
                       unlist(tissueDictionary$fetal[[1]]),
                       unlist(tissueDictionary$`stem-like`[[1]]),
                       unlist(tissueDictionary$`Neural progenitor/stem cells`[[1]]),
                       unlist(tissueDictionary$hESCs[[1]]),
                       unlist(tissueDictionary$iPSCs[[1]]),
                       unlist(tissueDictionary$`Other stem-cells`[[1]]),
                       unlist(tissueDictionary$`Other prenatal tissues`[[1]]),
                       unlist(tissueDictionary$fibroblasts[[1]]))
  lines <- readLines("Data/bulkRNASeq/cellosaurus.txt")
  
  # Get list of all cancer cell lines
  cancerCellVec <- c()
  for (i in 1:length(lines)) {
    line <- lines[i]
    grepRes1 <- grep(line, pattern = "^ID")
    grepRes2 <- grep(line, pattern = "^SY")
    grepRes3 <- grep(line, pattern = "^CA")
    if (length(grepRes1)) {
      tumors <- FALSE
      cellsNow <- gsub(x = line, pattern = "^ID\\s+",
                       replacement = "", perl = TRUE)
    }
    if (length(grepRes2)) {
      cellsNow2 <- gsub(x = line, pattern = "^SY\\s+",
                        replacement = "", perl = TRUE)
      cellsNow <- c(cellsNow, unlist(strsplit(cellsNow2, split = "; ")))
    }
    if (length(grepRes3)) {
      if (line == "CA   Cancer cell line") {
        cancerCellVec <- c(cancerCellVec, cellsNow)
      }
    }
  }
  
  toCheck <- colData$tissue[1]
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\+", replacement = "\\\\+")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\[", replacement = "\\\\[")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\]", replacement = "\\\\]")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\(", replacement = "\\\\(")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\)", replacement = "\\\\)")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\.", replacement = "\\\\.")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\^", replacement = "\\\\^")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\$", replacement = "\\\\$")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\?", replacement = "\\\\?")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\|", replacement = "\\\\|")
  cancerCellVec <- gsub(cancerCellVec, pattern = "\\*", replacement = "\\\\*")
  cancerCellVec2 <- gsub(cancerCellVec, pattern = "\\#", replacement = "")
  
  tissueUnique <- unique(colData$tissue)
  cancerCellVec1 <- sapply(cancerCellVec, FUN = function(x) {
    paste0(" ", x, " ")
  })
  cancerCellVec2 <- sapply(cancerCellVec, FUN = function(x) {
    paste0("^", x, " ")
  })
  cancerCellVec3 <- sapply(cancerCellVec, FUN = function(x) {
    paste0(" ", x, "$")
  })
  cancerCellVec4 <- sapply(cancerCellVec, FUN = function(x) {
    paste0("^", x, "$")
  })
  allGrep <- unique(c(cancerCellVec1, cancerCellVec2, cancerCellVec3, cancerCellVec4, cancerString))
  cancerCellList <- c()
  calculateGrep <- function(wordVec, patternVec) {
    # wordVec <-  tissueUnique[1:100]
    # patternVec <- cancerCellVec1[30344:3344]
    return(sapply(wordVec, FUN = grep, ignore.case = TRUE,
                  pattern = paste0(as.character(patternVec), collapse = "|"), perl = TRUE))
  }
  blocks <- 1000
  chunkSize <- ceiling(length(allGrep)/blocks)
  cancerCellList <- split(allGrep, ceiling(seq_along(allGrep)/chunkSize))
  registerDoParallel(cores=maxCores)
  resForEach <- foreach(k=1:blocks) %dopar% calculateGrep(tissueUnique, cancerCellList[[k]])
  cancerSamps <- unlist(lapply(resForEach, FUN = function(x) {
    return(names(x)[which(x != 0)])
  }))
  cancerSamps2 <- cancerSammps[grep(cancerSamps, invert = TRUE, ignore.case = TRUE, perl = TRUE,
                                    pattern = paste0(nonCancerString, collapse = "|"))]
  colData$disease <- "Normal"
  colData$disease[colData$tissue %in% cancerSamps2] <- "Tumor"
  load("Data/bulkRNASeq/poorlyDescribedTumorSamples.rda")
  colData$disease[colData$samples %in% missingTumorSamps] <- "Tumor"
  
  # Categorize remaining samples
  tissueDictNow <- tissueDict[c(1:31)]
  colData <- categorizeMetaData(metadata = colData,
                                cols = "tissue",
                                dictionary = tissueDictNow)
  rSums <- rowSums(colData[,c(36:39)])
  colDataMain <- unique(colData[rSums == 1,])
  colDataMain$Tissue <- "None"
  possibles <- colnames(colDataMain[,c(36:39)])
  resVec <- c()
  for (i in 1:length(colDataMain$samples)) {
    resVec <- c(resVec, possibles[which(colDataMain[i,c(36:39)] == 1)])
  }
  colDataMain$Tissue <- resVec
  
  # Get supplementary assignments
  rSums <- rowSums(colData[,c(9:35)])
  colDataSupp <- unique(colData[rSums == 1,])
  colDataSupp$Tissue <- "None"
  possibles <- colnames(colDataSupp[,c(9:35)])
  resVec <- c()
  for (i in 1:length(colDataSupp$samples)) {
    resVec <- c(resVec, possibles[which(colDataSupp[i,c(9:35)] == 1)])
  }
  resVec <- gsub(resVec, pattern = "repiratory", replacement = "respiratory")
  colDataSupp$Tissue <- resVec
  
  # Combine data
  colDataFinal <-  data.frame(
    "samples" = c(colDataMain$samples, colDataSupp$samples),
    "series" = c(colDataMain$series, colDataSupp$series),
    "tissue" = c(colDataMain$tissue, colDataSupp$tissue),
    "readsAligned" = c(colDataMain$reads_aligned, colDataSupp$reads_aligned),
    "disease" = c(colDataMain$disease, colDataSupp$disease),
    "TissueType" = c(colDataMain$Tissue, colDataSupp$Tissue), stringsAsFactors = F
  )
  colDataFinal <- colDataFinal[! duplicated(colDataFinal$sample),]
  
  # Final cleanup of tumor terms
  otherTumorTerms <- c("BeWo", "REH", "[a-zA-Z]oma$", "[a-zA-Z]oma ", "cancer", "tumor")
  colDataFinal$disease[grep(colDataFinal$tissue, ignore.case = TRUE, perl = TRUE,
                            pattern = paste0(otherTumorTerms, collapse = "|"))] <- "Tumor"
  
  # Remove Non-Ewing tumor samps
  tableDF <- as.data.frame(table(colDataFinal$TissueType, colDataFinal$disease))
  colDataFinal <- colDataFinal[colDataFinal$disease != "Tumor" |
                                 colDataFinal$TissueType == "ewing sarcoma",]
  # ES2 is not actually a EWS tumor...
  colDataFinal <- colDataFinal[grep(colDataFinal$tissue, ignore.case = FALSE, perl = TRUE,
                                    pattern = "ES2", invert = TRUE),]
  
  # Filter for read alignments > 10 million reads
  colDataFinal <- colDataFinal[colDataFinal$readsAligned > 10E6,]
  tableDF <- as.data.frame(table(colDataFinal$TissueType, colDataFinal$disease))
  
  save(colDataFinal, file = "Data/bulkRNASeq/colDataFinal.rda")
} else {
  cat("\ncolData File Found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/colDataFinal.rda")
}
cat("\n", timestamp2(), " Loading & filtering expression data...\n", sep = "")
if (! file.exists("Data/bulkRNASeq/fullRawCountsFiltered.rda")) {
  # Load and filter expression data
  expression <- h5read(dataFile, "data/expression")
  samples <- h5read(dataFile, "meta/Sample_geo_accession")
  genes <- h5read(dataFile, name = "meta/genes")
  rownames(expression) <- genes
  colnames(expression) <- samples
  colDataFinal <- colDataFinal[colDataFinal$samples %in% colnames(expression),]
  colDataFinalMain <- colDataFinal[colDataFinal$disease == "Normal" |
                                     colDataFinal$TissueType == "ewing sarcoma",]
  colDataFinalMain$TissueType[grep(colDataFinalMain$TissueType, ignore.case = T,
                                   pattern = "MSC|HSC|neural crest")] <- "stem-like"
  tableDF <- as.data.frame(table(colDataFinalMain$TissueType))
  colDataNow <- colDataFinalMain[colDataFinalMain$TissueType %in% tableDF$Var1[tableDF$Freq > 100],]
  expression <- expression[, which(colnames(expression) %in% colDataNow$samples)]
  colDataNow <- colDataNow[order(match(colDataNow$samples, colnames(expression))),]
  if (! all(colDataNow$samples == colnames(expression))) {
    stop("ColData samples are not identical to colnames of expression data...",
         " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
  }
  # Keep genes expressed in 10%+ of samples
  nonZeroCount <- apply(expression, 1, nonZeroSamps)
  expression <- expression[which(nonZeroCount > (length(colnames(expression)) * .10)),] 
  timestamp()
  cat("\nDone. Saving expression data...\n")
  save(expression, file = "Data/bulkRNASeq/fullRawCountsFiltered.rda")
} else if (file.exists("Data/bulkRNASeq/fullVSTCounts.rda")) {
  cat("\nVST-transformed expression data found... Skipping this step...\n")
} else {
  cat("\nFiltered expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
}

cat("\n", timestamp2(), " Applying Variance Stabilizing Transform to dataset...\n", sep = "")
UMAP <- FALSE
if (! file.exists("Data/bulkRNASeq/fullVSTCounts.rda")) {
  # Make data Homoscedastic with VST 
  vsd <- vst(expression, nsub = 1000)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/bulkRNASeq/fullVSTCounts.rda")
} else if (file.exists("Data/bulkRNASeq/fullUMAPData.rda")) {
  UMAP <- TRUE
  cat("\nUMAP results found... Skipping this step...\n")
} else {
  cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}

if (! file.exists("Data/bulkRNASeq/fullUMAPData.rda") & ! UMAP) {
  cat("\n", timestamp2(), " Calculating UMAP...\n", sep = "")
  # Calculate highly-variable genes (HVGs)
  rv <- matrixStats::rowVars(vsd)
  names(rv) <- rownames(vsd)
  rv <- rv[order(rv, decreasing = T)]
  hvgs <- rv[c(1:10000)]
  vsdHVG <- vsd[names(hvgs),]
  rm(vsd)
  gc()
  colDataNow <- colDataFinal[colDataFinal$samples %in% colnames(vsdHVG),]
  colDataNow <- colDataNow[order(match(colDataNow$samples, colnames(vsdHVG))),]
  if (! all(colDataNow$samples == colnames(vsdHVG))) {
    stop("ColData samples are not identical to colnames of expression data...",
         " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
  }
  pltData <- colDataNow
  
  # Calculate PCA
  if (! file.exists( "Data/bulkRNASeq/fullPCData.rda")) {
    pcDataFull <- prcomp(t(vsdHVG))
    pcData <- pcDataFull$x
    pcData <- pcData[,c(1:100)]
    save(pcData, file = "Data/bulkRNASeq/fullPCData.rda")
  } else {
    load("Data/bulkRNASeq/fullPCData.rda")
  }
  pltData$PC1 <- pcData[,c(1)]
  pltData$PC2 <- pcData[,c(2)]
  # Calculate Neighbors
  neighbors <- FindNeighbors(pcData)
  # Calculate clustering
  clusters <- FindClusters(neighbors$snn, resolution = .015)
  pltData$cluster <- clusters[,c(1)]
  # Calculate UMAP embedding  
  nnNow <- floor(length(colnames(vsdHVG)) * .2)
  umapData <- umap(t(vsdHVG), verbose = T,
                   n_neighbors = nnNow, 
                   min_dist = 1, pca = 100)
  pltData$UMAP_1 <- umapData[,c(1)]
  pltData$UMAP_2 <- umapData[,c(2)]
  save(pltData, file = "Data/bulkRNASeq/fullUMAPData.rda")
} else {
  if (UMAP) {
    cat("\nUMAP data found! Loading and continuing to next step...\n")
  }
  load("Data/bulkRNASeq/fullUMAPData.rda")
}

# Find EWS cell type markers
cat("\n", timestamp2(), " Calculating cell type marker genes for EWS...\n", sep = "")
if (! file.exists("Data/bulkRNASeq/ewsMarkerGenes.rda")) {
  # Load full UMAP/VST data
  load("Data/bulkRNASeq/fullUMAPData.rda")
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  all(pltData$samples == colnames(expression))
  srt <- CreateSeuratObject(expression, meta.data = pltData)
  srt <- NormalizeData(srt)
  srt <- ScaleData(srt, block.size = 10000)
  srt$TissueType <- pltData$TissueType
  Idents(srt) <- srt$TissueType
  srtMarkers <- FindMarkers(srt, ident.1 = "ewing sarcoma", logfc.threshold = 0,
                            only.pos = TRUE, random.seed = 42)
  ewsMarkerGenes <- srtMarkers
  save(ewsMarkerGenes, file = "Data/bulkRNASeq/ewsMarkerGenes.rda")
  # save(srt, file = "Data/bulkRNASeq/srt.rda")
} else {
  cat("\nCell type marker genes found! Loading and continuing...\n")
  load("Data/bulkRNASeq/ewsMarkerGenes.rda")
}

# Plot bulk RNA-Seq figure 1
cat("\n", timestamp2(), " Plotting first bulk RNA-Seq figures...\n", sep = "")
if (! file.exists("Data/bulkRNASeq/pltDataUMAP_Fig1.rda")) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
  load("Data/bulkRNASeq/fullUMAPData.rda")
  # Classify tissue types
  tissueDictNow <- jsonlite::read_json("Data/bulkRNASeq/tissueDictionary.json")
  pltDataToClass <- unique(pltData[,c(1:3)])
  pltDataToClass <- categorizeMetaData(pltDataToClass,
                                       cols = c("tissue", "TissueType"), 
                                       dictionary = tissueDictNow[c(1:6,8:15,17:31, 33:39)])
  rSums <- rowSums(pltDataToClass[,c(4:39)])
  toClass <- pltDataToClass[rSums == 1,]
  classes <- c()
  toClass$Tissue <- "None"
  possibles <- colnames(toClass[,c(4:39)])
  resVec <- c()
  for (i in 1:length(toClass$samples)) {
    resVec <- c(resVec, possibles[which(toClass[i,c(4:39)] == 1)])
  }
  toClass$Tissue <- resVec
  table(toClass$Tissue)
  classMap <- toClass[,c(1, 2, 3, 40)]
  pltData <- pltData[,c(-13)]
  catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
  catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
  catClass$Tissue[is.na(catClass$Tissue) & 
                    catClass$TissueType == "stem-like"] <- "Other stem-cells"
  catClass$Tissue[is.na(catClass$Tissue) & 
                    catClass$TissueType == "MSCs"] <- "MSCs"
  catClass$Tissue[is.na(catClass$Tissue) & 
                    catClass$TissueType == "HSCs"] <- "HSCs"
  table(catClass$Tissue)
  pltData$TissueFinal <- catClass$Tissue
  pltData <- pltData[pltData$TissueFinal != "Non-Ewing tumor",]
  pltData$TissueFinal[is.na(pltData$TissueFinal)] <- "Mixed/other"
  vsd <- vsd[, colnames(vsd) %in% pltData$samples]
  pltData <- pltData[pltData$samples %in% colnames(vsd),]
  all(pltData$samples == colnames(vsd))
  dim(vsd)
  pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                   "iPSCs", "hESCs")] <- 
    str_to_sentence(pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                                     "iPSCs", "hESCs")])
  # Get EWS score
  load("Data/bulkRNASeq/ewsMarkerGenes.rda")
  ewsMarkerGenes$p_val_adj[ewsMarkerGenes$p_val_adj == 0] <- .Machine$double.xmin
  pltData$Cluster <- pltData$cluster
  ewingMarkers <- rownames(ewsMarkerGenes)[ewsMarkerGenes$avg_logFC > .58 &
                                             ewsMarkerGenes$p_val_adj < 1E-50]
  ewingMarkers <- ewingMarkers[ewingMarkers %in% rownames(vsd)]
  ewsMarkerGenes$p_val_adj[ewsMarkerGenes$p_val_adj == 0] <- .Machine$double.xmin
  ewsMarkerGenes$pAdj <- -log10(ewsMarkerGenes$p_val_adj)
  ewsMarkerGenes$geneName <- rownames(ewsMarkerGenes)
  ewsMarkerGenes$Marker <- FALSE
  ewsMarkerGenes$Marker[ewsMarkerGenes$avg_logFC > .58 &
                          ewsMarkerGenes$p_val_adj < 1E-50] <- TRUE
  doMarksDF <- ewsMarkerGenes[ewsMarkerGenes$avg_logFC > 1.5 &
                                ewsMarkerGenes$p_val_adj < 1E-50,]
  g1 <- ggplot(ewsMarkerGenes, mapping = aes_string(x = "avg_logFC", y = "pAdj",
                                                    color = "Marker")) +
    geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
    rremove("legend") + 
    ylab("P adj value (-log10)") +
    xlab("Log2 fold change") +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,
                                                   (max(ewsMarkerGenes$pAdj)*1.05))) +
    scale_x_continuous(expand = c(0,0), limits = c(0, (max(ewsMarkerGenes$avg_logFC)*1.05))) +
    scale_color_manual(values=c( "#C4BFBF", "#CC0909")) +
    labs(color = "EWS marker") +
    geom_text_repel(data = doMarksDF,segment.colour = "black", 
                    color = "black", seed = 42, #nudge_x = .15,
                    # force = 2, point.padding = .5,
                    mapping = aes_string(x = "avg_logFC", y = "pAdj",
                                         label = "geneName"))
  ggsave(g1, filename = paste0("Figures_v2/FigS1A.png"), height = 7.5, 
         width = 8)
  g1 <- ggplot(ewsMarkerGenes, mapping = aes_string(x = "avg_logFC", y = "pAdj",
                                                    color = "Marker")) +
    geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
    ylab("P adj value (-log10)") +
    xlab("Log2 fold change") +
    theme(legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,
                                                   (max(ewsMarkerGenes$pAdj)*1.05))) +
    scale_x_continuous(expand = c(0,0), limits = c(0, (max(ewsMarkerGenes$avg_logFC)*1.05))) +
    scale_color_manual(values=c( "#C4BFBF", "#CC0909")) +
    labs(color = "EWS marker") +
    geom_text_repel(data = doMarksDF,segment.colour = "black", 
                    color = "black", seed = 42, #nudge_x = .15,
                    # force = 2, point.padding = .5,
                    mapping = aes_string(x = "avg_logFC", y = "pAdj",
                                         label = "geneName"))
  ggsave(g1, filename = paste0("Figures_v2/FigS1A_legend.png"), height = 7.5, 
         width = 10.64)
  
  # Marker genes analysis
  marksNow <- ewsMarkerGenes$geneName
  marksNowTERM <- marksNow[which(marksNow %in% TERM2GENE_EWS$gene_symbol)]
  marksNowBIND <- marksNow[which(marksNow %in% TERM2GENE_EWS$gene_symbol[grep(x = TERM2GENE_EWS$gs_name, pattern = "TARGET|ACTIVATE|BOUND")])]
  marksNowBIND <- marksNowBIND[which(marksNowBIND %in% TERM2GENE_EWS$gene_symbol[grep(x = TERM2GENE_EWS$gs_name, pattern = "AYNAUD|RIGGI|SILIGAN")])]
  # Direct EWSR1-FLI1 targets
  compList <- list(
    "EWS Markers" = ewsMarkerGenesFinal$geneName,
    "Bound by EWSR1-FLI1" = marksNowBIND,
    "Ewing sarcoma gene sets" = marksNowTERM
  )
  calculate.overlap.and.pvalue(compList[[1]], compList[[2]], lower.tail = F,
                               total.size = length(intersect(TERM2GENE$gene_symbol, 
                                                      ewsMarkerGenes$geneName)))
  gd <- venn.diagram(compList[c(1,2)], filename = NULL,cat.dist= c(.25, .25),
                     fill = c("Firebrick", "skyblue"))
  plot.new()
  dev.off()
  grid.draw(gd)
  gd <- grid.grab()
  plt <- as.ggplot(gd)
  dev.off()
  ggsave(plt, filename = "Figures_v2/Venn_Diagram_BoundByEF1_EWSMarkers.pdf", height = 4, width = 6)
  calculate.overlap.and.pvalue(compList[[1]], compList[[3]], lower.tail = F,
                               total.size = length(intersect(TERM2GENE$gene_symbol, 
                                                             ewsMarkerGenes$geneName)))
  gd <- venn.diagram(compList[c(1,3)], filename = NULL,cat.dist= c(.25, .25),
                     fill = c("Firebrick", "palegoldenrod"))
  plot.new()
  dev.off()
  grid.draw(gd)
  gd <- grid.grab()
  plt <- as.ggplot(gd)
  dev.off()
  ggsave(plt, filename = "Figures_v2/Venn_Diagram_EwingGeneSets_EWSMarkers.pdf", height = 4, width = 6)
  
  
  ewsMarkerGenes2 <- ewsMarkerGenesFinal
  ewsMarkerGenes2$Group <- "Other"
  ewsMarkerGenes2$Group[ewsMarkerGenes2$geneName %in% ewsMarkerGenesFinal] <- "EWS Marker"
  ewsMarkerGenes2$Group[ewsMarkerGenes2$geneName %in% marksNowBIND] <- "EWSR1-FLI1 Target"
  ewsMarkerGenes2 <- ewsMarkerGenes2[order(ewsMarkerGenes2$Group),]
  
  g1 <- ggplot(ewsMarkerGenes2, mapping = aes_string(x = "avg_logFC", y = "pAdj",
                                                    color = "Group")) +
    geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
    ylab("P adj value (-log10)") +
    xlab("Log2 fold change") +
    theme(legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    scale_y_continuous(expand = c(0,0), limits = c(0,
                                                   (max(ewsMarkerGenes2$pAdj)*1.05))) +
    scale_x_continuous(expand = c(0,0), limits = c(0, (max(ewsMarkerGenes2$avg_logFC)*1.05))) +
    scale_color_manual(values=rev(c( "#C4BFBF", "#00B0F6", "#CC0909")) )
  ggsave(g1, filename = paste0("Figures_v2/EWS_Markers_and_Binding_Targets.png"), height = 7.5, 
         width = 10.64)
  
  
  
  # TERM2GENE <- getTERM2GENE("complex")
  # ranks <- -log10(ewsMarkerGenes$p_val_adj) * ewsMarkerGenes$avg_logFC
  # names(ranks) <- rownames(ewsMarkerGenes)
  # GSEARes <- myGSEA(ranks, TERM2GENE = TERM2GENE)
  # eres <- GSEARes$eres
  
  pltData$ewsScore <- colMedians(vsd[ewingMarkers,])
  pltData$Cluster <- pltData$cluster
  pltData$Tissue <- pltData$TissueFinal
  save(pltData, file = "Data/bulkRNASeq/pltDataUMAP_Fig1.rda")
} else {
  load("Data/bulkRNASeq/pltDataUMAP_Fig1.rda")
}

# Get Ewing marker gene table
if (! "vsd" %in% ls()) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}
ewsMarkerGenes$p_val_adj[ewsMarkerGenes$p_val_adj == 0] <- .Machine$double.xmin
ewingMarkers <- rownames(ewsMarkerGenes)[ewsMarkerGenes$avg_logFC > .58 &
                                           ewsMarkerGenes$p_val_adj < 1E-50]
ewingMarkers <- ewingMarkers[ewingMarkers %in% rownames(vsd)]
ewsMarkerGenes$pAdj <- -log10(ewsMarkerGenes$p_val_adj)
ewsMarkerGenes$geneName <- rownames(ewsMarkerGenes)
ewsMarkerGenes$Marker <- FALSE
ewsMarkerGenes$Marker[ewsMarkerGenes$avg_logFC > .58 &
                        ewsMarkerGenes$p_val_adj < 1E-50] <- TRUE
ewsMarkerGenes <- ewsMarkerGenes[ewsMarkerGenes$Marker,]
ewsMarkerGenesFinal <- ewsMarkerGenes[,c(7, 2, 1, 5)]
ewsMarkerGenesFinal <- ewsMarkerGenesFinal[order(ewsMarkerGenesFinal$avg_logFC, decreasing = TRUE),]
write.table(ewsMarkerGenesFinal, file = "Tables/ewsMarkerGenes.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Write colData table
colDataFinal <- pltData[,c(-6, -12)]
colDataFinal <- colDataFinal[, c(2, 1, 3, 4, 13, 5, 12, 11, 6, 7, 9, 10)]
colnames(colDataFinal) <- c("Series", "Sample", "originalTissue", "readsAligned", "Tissue",
                            "Disease", "ewsScore", "Cluster",
                            "PC1", "PC2", "UMAP_1", "UMAP_2")
write.table(colDataFinal, file = "Tables/colDataBulk.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Clear old objects
suppressWarnings(rm("vsdHVG"))
suppressWarnings(rm("srt"))
suppressWarnings(rm("pcDataFull"))
suppressWarnings(rm("resForEach"))
invisible(gc())

# Make color map
gg_color <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gCols <- gg_color(4)
colMapCategory <- c("Ewing sarcoma" = gCols[1],
                    "Mesodermal" = gCols[2],
                    "Neuroectodermal" = gCols[3],
                    "Pluripotent" = gCols[4],
                    "Mixed/other" = "#C1C1C1")
gCols <- gg_color(7)
colMapCluster <- c("2" = gCols[1],
                   "4" = gCols[2],
                   "28" = gCols[3],
                   "17" = gCols[4],
                   "19" = gCols[5],
                   "16" = gCols[6],
                   "24" = gCols[7])
gCols <- gg_color(14)
colMapTissueNow <- c("Ewing sarcoma" = gCols[1],
                     "Brain" = gCols[2],
                     "Cardiac" = gCols[3],
                     "Endothelial" = gCols[4],
                     "Fibroblasts" = gCols[5],
                     "hESCs" = gCols[6],
                     "HSCs" = gCols[7],
                     "iPSCs" = gCols[8],
                     "MSCs" = gCols[9],
                     "Muscle" = gCols[10],
                     "Neural crest cells" = gCols[11],
                     "Neural progenitor/stem cells" = gCols[12],
                     "Other prenatal tissues" = gCols[13],
                     "Other stem-cells" = gCols[14],
                     "Mixed/other" = "#C1C1C1")

# Plot clusters and tissue types
DimPlot2(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                              color = "Cluster"), 
         plotName = "Fig1A")
DimPlot2(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                              color = "Tissue"), 
         plotName = "Fig1B")
# EWS Markers highlight EWS cells -- some other tissues have higher expression of these
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "ewsScore")) +
  geom_point(size = .4) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + 
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "OrRd"))(255))
ggsave(g1, filename = paste0("Figures_v2/Fig1C.png"), height = 7.5, 
       width = 8)
g1 <- ggplot(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                  color = "ewsScore")) +
  geom_point(size = .4) + theme_pubr(border = T, base_size = 22) +
  # guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  labs(color = "EWS marker\nexpression") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "OrRd"))(255))
ggsave(g1, filename = paste0("Figures_v2/Fig1C_legend.png"), height = 7.5, 
       width = round((8 * 1.33), digits = 2))


# What clusters and tissue types show the highest EWS score?
fillVar <- "TissueFinal"
labDF <- pltData
labDF$cluster <- paste0("Cluster_", labDF$cluster)
orderDF1 <- labDF %>% dplyr::group_by(eval(parse(text = fillVar))) %>%
  select(.data$ewsScore) %>% summarise_all(median)
colnames(orderDF1)[1] <- fillVar
orderDF1 <- orderDF1[order(orderDF1$ewsScore, decreasing = T),]
labDF <- labDF[order(match(labDF[,fillVar, drop = T], orderDF1[,fillVar, drop = T])),]
levels <- rev(unique(labDF[,fillVar]))
labDF[,fillVar] <- factor(labDF[,fillVar], levels =levels)
pal <- colorRampPalette(brewer.pal(9, "OrRd"))(length(unique(labDF[,c(fillVar)])))
g1 <- ggplot(labDF, mapping = aes_string(x = fillVar, y = "ewsScore", fill = fillVar)) +
  geom_boxplot() +  
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme_pubr(border = T, base_size = 22) +
  ylab(paste0("EWS marker expression")) +
  rremove("ylab") +
  rremove("legend") + 
  scale_fill_manual(values = pal) +
  # scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "OrRd"))(255)) +
  ggpubr::rotate()
ggsave(g1, filename = paste0("Figures_v2/Fig1D.png"), height = 10, width = 10)

# Clusters
fillVar <- "cluster"
labDF <- pltData
labDF$cluster <- paste0("Cluster_", labDF$cluster)
orderDF2 <- labDF %>% dplyr::group_by(eval(parse(text = fillVar))) %>%
  select(.data$ewsScore) %>% summarise_all(median)
colnames(orderDF2)[1] <- fillVar
orderDF2 <- orderDF2[order(orderDF2$ewsScore, decreasing = T),]
labDF <- labDF[order(match(labDF[,fillVar, drop = T], orderDF2[,fillVar, drop = T])),]
levels <- rev(unique(labDF[,fillVar]))
labDF[,fillVar] <- factor(labDF[,fillVar], levels =levels)
pal <- colorRampPalette(brewer.pal(9, "OrRd"))(length(unique(labDF[,c(fillVar)])))
g1 <- ggplot(labDF, mapping = aes_string(x = fillVar, y = "ewsScore", fill = fillVar)) +
  geom_boxplot() +  
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme_pubr(border = T, base_size = 22) +
  ylab(paste0("EWS marker expression")) +
  rremove("ylab") +
  rremove("legend") + 
  scale_fill_manual(values = pal) +
  # scale_fill_gradientn(colors = colorRampPalette(brewer.pal(9, "OrRd"))(255)) +
  ggpubr::rotate()
ggsave(g1, filename = paste0("Figures_v2/Fig1E.png"), height = 10, width = 8)

# Heatmap showing tissues in clusters
labDF <- pltData
# labDF$cluster <- paste0("Cluster_", labDF$cluster)
tableDF <- as.data.frame(table(labDF$TissueFinal, labDF$cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Blues")))(paletteLength)
resMat <- resMat[orderDF2$cluster,orderDF1$TissueFinal]
annRow <- data.frame(row.names = orderDF2$cluster, stringsAsFactors = F,
                     "EWS_score" = orderDF2$ewsScore)
annCol <- data.frame(row.names = orderDF1$TissueFinal, stringsAsFactors = F,
                     "EWS_score" = orderDF1$ewsScore)
annColor <- colorRampPalette(brewer.pal(9, "OrRd"))(255)
ph <- pheatmap(resMat, cluster_rows = F, cluster_cols = F, 
               annotation_col = annCol,
               annotation_colors = list("EWS_score" = annColor),
               annotation_names_row= F, silent = TRUE,
               annotation_row = annRow, annotation_names_col = F,
               labels_col = as.character(colnames(resMat)),
               color = heatColor, angle_col = 315, 
               margins = c(5, 19), fontsize = 13.5)
g1 <- as.ggplot(ph[[4]])
ggsave(g1, filename = paste0("Figures_v2/FigS1B.png"), 
       height = 8, width = 12)

# Calculate PHATE on bulk RNA-Seq data 
cat("\n", timestamp2(), " Calculating and tuning PHATE with plotting...\n", sep = "")
if (! file.exists("Data/bulkRNASeq/pltDataNonEWS_PHATEFinal.rda")) {
  # Choosing the top clusters
  topSamps <- pltData$samples[pltData$cluster %in% c(as.numeric(names(colMapCluster)))]
  pltDataTop <- pltData[pltData$samples %in% topSamps,]
  vsdTop <- vsd[,topSamps]
  # Calculate PHATE
  kNow <- 350
  phateNow <- phate(t(vsdTop), knn = kNow, 
                    n.jobs = -1, seed = 42, ndim = 3)
  pltDataTop$PHATE_1 <- phateNow$embedding[,c(1)]
  pltDataTop$PHATE_2 <- phateNow$embedding[,c(2)]
  pltDataTop$PHATE_3 <- phateNow$embedding[,c(3)]
  pltDataTop$Cluster <- pltDataTop$cluster
  DimPlot2(data = pltDataTop, mapping = aes_string(x = "PHATE_2", y = "PHATE_3",
                                                   color = "Cluster"), 
           plotName = "FigS1G")
  rmSamps <- pltDataTop$samples[pltDataTop$PHATE_2 > 0.007 |
                                  (pltDataTop$PHATE_1 < 0 & pltDataTop$PHATE_2 > 0)]
  pltDataTop$Remove <- FALSE
  pltDataTop$Remove[pltDataTop$samples %in% rmSamps] <- TRUE
  DimPlot2(data = pltDataTop, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                   color = "Remove"), 
           plotName = "FigS1H")
  
  # Next filter
  pltDataTop2 <- pltDataTop[! pltDataTop$samples %in% rmSamps,]
  vsdTop2 <- vsdTop[, colnames(vsdTop) %in% pltDataTop2$samples]
  phateNow <- phate(t(vsdTop2), knn = kNow, 
                    n.jobs = -1, seed = 42, ndim = 3)
  pltDataTop2$PHATE_1 <- phateNow$embedding[,c(1)]
  pltDataTop2$PHATE_2 <- phateNow$embedding[,c(2)]
  pltDataTop2$PHATE_3 <- phateNow$embedding[,c(3)]
  pltDataTop2$Cluster <- pltDataTop2$cluster
  DimPlot2(data = pltDataTop2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                    color = "Cluster"), 
           plotName = "FigS1I")
  
  rmSamps <- pltDataTop2$samples[pltDataTop2$PHATE_2 > 0.0025 & 
                                   pltDataTop2$PHATE_1 > 0 ]
  pltDataTop2$Remove <- FALSE
  pltDataTop2$Remove[pltDataTop2$samples %in% rmSamps] <- TRUE
  DimPlot2(data = pltDataTop2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                    color = "Remove"), 
           plotName = "FigS1J")
  # Last filter
  pltDataTopFinal <- pltDataTop2[! pltDataTop2$samples %in% rmSamps,]
  vsdTopFinal <- vsdTop2[, colnames(vsdTop2) %in% pltDataTopFinal$samples]
  
  # Check PHATE QC one more time
  kNow <- 350
  phateNow <- phate(t(vsdTopFinal), knn = kNow, 
                    n.jobs = -1, seed = 42, ndim = 3)
  pltDataTopFinal$PHATE_1 <- phateNow$embedding[,c(1)]
  pltDataTopFinal$PHATE_2 <- phateNow$embedding[,c(2)]
  pltDataTopFinal$PHATE_3 <- phateNow$embedding[,c(3)]
  pltDataTopFinal$Cluster <- pltDataTopFinal$cluster
  DimPlot2(data = pltDataTopFinal, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                        color = "Cluster"), 
           plotName = "FigS1K")
  
  ## Main PHATE analysis
  # Categorize germ layer
  pltDataNow <- pltDataTopFinal
  pltDataNow$TissueNow <- "Mixed/other"
  intTissue <- c("Ewing sarcoma", "Fibroblasts", 
                 "Brain", "hESCs", "iPSCs", "Neural progenitor/stem cells",
                 "HSCs", "Endothelial", "Other prenatal tissues",
                 "MSCs", "Neural crest cells", "Other stem-cells")
  pltDataNow$TissueNow[pltDataNow$TissueFinal %in% intTissue] <- 
    pltDataNow$TissueFinal[pltDataNow$TissueFinal %in% intTissue]
  MesodermalTissues <- c("Adipose", "Bone", "Fibroblasts", 
                          "MSCs", "Endothelial", "Cardiac",
                          "Muscle", "Male reproductive")
  ectodermalTissues <- c("Brain", "Retina",
                         "Neural progenitor/stem cells")
  PluripotentTissues <- c("hESCs", "iPSCs", 
                         "Other prenatal tissues", 
                         "Neural crest cells")
  pltDataNow$Category <- "Mixed/other"
  pltDataNow$Category[pltDataNow$TissueFinal %in%
                        MesodermalTissues] <- "Mesodermal"
  pltDataNow$Category[pltDataNow$TissueFinal %in%
                        ectodermalTissues] <- "Neuroectodermal"
  pltDataNow$Category[pltDataNow$TissueFinal %in%
                        PluripotentTissues] <- "Pluripotent"
  pltDataNow$Category[pltDataNow$TissueFinal %in%
                        "Ewing sarcoma"] <- "Ewing sarcoma"
  
  # Run PHATE
  nnList <- seq(50, 1000, by = 50)
  dir.create("Figures_v2/PHATEoptimization", showWarnings = FALSE)
  for (i in 1:length(nnList)) {
    kNow <- nnList[i]
    print(kNow)
    if (! file.exists(paste0("PHATEoptimization/kNN_", kNow, ".png"))) {
      phateNow <- phate(t(vsdTopFinal), knn = kNow,
                        n.jobs = -1, seed = 42, ndim = 2)
      pltDataNow$PHATE_1 <- phateNow$embedding[,c(1)]
      pltDataNow$PHATE_2 <- phateNow$embedding[,c(2)]
      pltDataNow$Cluster <- pltDataNow$cluster
      DimPlot2(data = pltDataNow, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                       color = "Cluster"),
               colorMap = colMapCluster,
               plotName = paste0("PHATEoptimization/kNN_", kNow))
    } else {
      print("Already calculated for this value.")
    }
  }
  
  kNow <- 750
  phateNow <- phate(t(vsdTopFinal), knn = kNow, 
                    n.jobs = -1, seed = 42, ndim = 2)
  pltDataNow$PHATE_1 <- phateNow$embedding[,c(1)]
  pltDataNow$PHATE_2 <- phateNow$embedding[,c(2)]
  pltDataNow$Cluster <- pltDataNow$cluster
  pltDataNow$Category <- factor(pltDataNow$Category, levels = names(colMapCategory))
  pltDataNow$TissueNow <- factor(pltDataNow$TissueNow, levels = names(colMapTissueNow))
  pltDataNonEWS <- pltDataNow[pltDataNow$TissueFinal != "Ewing sarcoma",]
  vsdTopFinalNonEWS <- vsdTopFinal[, colnames(vsdTopFinal) %in% pltDataNonEWS$samples]
  phateNonEWS <- phate(t(vsdTopFinalNonEWS), knn = kNow, 
                       n.jobs = -1, seed = 42, ndim = 2)
  pltDataNonEWS$PHATE_1 <- phateNonEWS$embedding[,c(1)]
  pltDataNonEWS$PHATE_2 <- phateNonEWS$embedding[,c(2)]
  pltDataNonEWS$Cluster <- pltDataNonEWS$cluster
  pltDataNonEWS$Category <- droplevels(pltDataNonEWS$Category)
  pltDataNonEWS$TissueNow <- droplevels(pltDataNonEWS$TissueNow)
  save(pltDataNow, file = "Data/bulkRNASeq/pltDataNow_PHATEFinal.rda")
  save(pltDataNonEWS, file = "Data/bulkRNASeq/pltDataNonEWS_PHATEFinal.rda")
} else {
  cat("\nPHATE data found! Loading and continuing...\n")
  load("Data/bulkRNASeq/pltDataNow_PHATEFinal.rda")
  load("Data/bulkRNASeq/pltDataNonEWS_PHATEFinal.rda")
}
pltDataNow$Category <- as.character(pltDataNow$Category)
pltDataNow$Category[pltDataNow$Category ==  "Mesenchymal"] <- "Mesodermal"
pltDataNow$Category[pltDataNow$Category ==  "Primordial"] <- "Pluripotent"
pltDataNow$Category[pltDataNow$Category ==  "Neuro-ectodermal"] <- "Neuroectodermal"
pltDataNow$Category[pltDataNow$TissueNow ==  "Neural crest cells"] <- "Neuroectodermal"
pltDataNow$Category <- factor(pltDataNow$Category, levels = c("Ewing sarcoma", "Mesodermal",
                                                              "Neuroectodermal", "Pluripotent", "Mixed/other"))
pltDataNonEWS$Category <- as.character(pltDataNonEWS$Category)
pltDataNonEWS$Category[pltDataNonEWS$Category ==  "Mesenchymal"] <- "Mesodermal"
pltDataNonEWS$Category[pltDataNonEWS$Category ==  "Primordial"] <- "Pluripotent"
pltDataNonEWS$Category[pltDataNonEWS$Category ==  "Neuro-ectodermal"] <- "Neuroectodermal"
pltDataNonEWS$Category[pltDataNonEWS$TissueNow ==  "Neural crest cells"] <- "Neuroectodermal"
pltDataNonEWS$Category <- factor(pltDataNonEWS$Category)
pltDataNonEWS$Category <- factor(pltDataNonEWS$Category, levels = c("Ewing sarcoma", "Mesodermal",
                                                              "Neuroectodermal", "Pluripotent", "Mixed/other"))

## PHATE visualizations
# Save PHATE embedding data to table
phateDataSave <- pltDataNow[,c(2, 1, 3, 20, 21, 13, 16, 17, 18)]
colnames(phateDataSave)[c(1:4)] <- c("Series", "Samples", "originalTissue", "Tissue")
write.table(phateDataSave, file = "Tables/phateDataSamples.tsv", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Show clusters w and w/out EWS
DimPlot2(data = pltDataNow, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                 color = "Cluster"), 
         colorMap = colMapCluster,
         plotName = "FigS1Li")
DimPlot2(data = pltDataNow, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                 color = "Category"), 
         colorMap = colMapCategory,
         plotName = "Fig1_PHATE_Categories")
pltDataNow2 <- pltDataNow
pltDataNow2$Category[pltDataNow2$Category != "Ewing sarcoma"] <- "Mixed/other"
pltDataNow2$Category <- droplevels(pltDataNow2$Category)
pltDataNow2 <- pltDataNow2[order(match(pltDataNow2$Category, levels(pltDataNow2$Category)), decreasing = TRUE),]
DimPlot2(data = pltDataNow2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                  color = "Category"), 
         colorMap = colMapTissueNow,
         plotName = "FigS1Lii")
pltDataNow2 <- pltDataNow
pltDataNow2$Category[pltDataNow2$Category == "Ewing sarcoma"] <- "Mixed/other"
pltDataNow2 <- pltDataNow2[order(match(pltDataNow2$Category, levels(pltDataNow2$Category)), decreasing = TRUE),]
DimPlot2(data = pltDataNow2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                  color = "Category"), 
         colorMap = colMapCategory,
         plotName = "FigS1Liii")
DimPlot2(data = pltDataNonEWS, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                    color = "Cluster"), 
         colorMap = colMapCluster,
         plotName = "FigS1M")
DimPlot2(data = pltDataNonEWS, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                    color = "Category"), 
         colorMap = colMapCategory,
         plotName = "FigS1M_category")

# Compare PHATE_1 for different EWS conditions
phateEWS <- pltDataNow[pltDataNow$TissueFinal == "Ewing sarcoma",]
EWSInfo <- read.csv("Data/bulkRNASeq/EWS_bulkRNASeq_sampleInfo.csv", 
                    stringsAsFactors = FALSE)
write.table(EWSInfo, file = "Tables/ewsSampleInfoDetails.tsv", 
            quote = FALSE, row.names = FALSE, sep = "\t")
phateEWS <- merge(x=phateEWS, all = TRUE,
                  y = EWSInfo, by = "samples")
phateEWS <- phateEWS[,c(22:24, 16:17)]
GSEToPlot <- unique(phateEWS$SeriesGSE)
pltList <- list()
for (i in 1:length(GSEToPlot)) {
  GSENow <- GSEToPlot[i]
  dataNow <- phateEWS[phateEWS$SeriesGSE == GSENow,]
  dataNow <- dataNow[! is.na(dataNow$PHATE_1),]
  dataNow$Condition <- paste0(dataNow$Cell, " ", dataNow$Condition)
  dataNow$PHATE_1Score <- range01(dataNow$PHATE_1)
  dataNow$Condition <- factor(dataNow$Condition, 
                              levels = rev(unique(dataNow$Condition)))
  pltNow <- ggplot(data = dataNow, aes_string(x = "Condition",
                                              y = "PHATE_1Score")) +
    geom_boxplot() +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") + ylab("Scaled PHATE_1 Score") +
    scale_fill_manual(values=colors) +
    labs(title = GSENow) + rremove('ylab') +
    theme(axis.text.y = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    scale_y_continuous(expand = c(0,0), limits = c(-.1,(max(dataNow$PHATE_1Score)*1.1))) +
    rremove("legend") + ggpubr::rotate()
  pltList[[i]] <- pltNow
}
ga <- ggarrange(plotlist = pltList, ncol = 5, nrow = 5, align = "hv")
ggsave(ga, height = 25, width = 40, filename = "Data/bulkRNASeq/PHATE_1_comps.png")

# Compare PHATE_1 for EWSFLI1 KD experiments
GSEToPlot <- c("GSE53066", "GSE60949", "GSE61950", 
               "GSE103843", "GSE92739",
               "GSE122537")
dataNow <- phateEWS[phateEWS$SeriesGSE %in% GSEToPlot,]
pltList <- list()
compListList <- list(
  list(c("A673 shEF1", "A673 shCTR")),
  list(c("A673 shEF1", "A673 shCTR")),
  list(c("SKNMC shFLI1", "SKNMC shCTR")),
  list(c("TC32 siFLI1", "TC32 siCTR")),
  list(c("A673 siCTR-shEF1", "A673 siCTR"),
       c("SKNMC shEF1", "SKNMC shCTR")),
  list(c("A673 shEF1","A673 shCTR"))
  
)
levelList <- list(
  c("A673 shCTR", "A673 shEF1"),
  c("A673 shCTR", "A673 shEF1", "A673 shEWSAT1"),
  c("SKNMC shCTR", "SKNMC shFLI1", "A673 shCTR", "A673 shFLI1"),
  c("TC32 siCTR", "TC32 siFLI1", "TC32 siEWSR1"),
  c("A673 siCTR", "A673 siCTR-shEF1", "SKNMC shCTR", "SKNMC shEF1"),
  c("A673 shCTR", "A673 shEF1", "A673 shEF1-oeEF1")
)
for (i in 1:length(GSEToPlot)) {
  GSENow <- GSEToPlot[i]
  dataNow <- phateEWS[phateEWS$SeriesGSE == GSENow,]
  dataNow <- dataNow[! is.na(dataNow$PHATE_1),]
  dataNow$Condition <- paste0(dataNow$Cell, " ", dataNow$Condition)
  dataNow <- dataNow[dataNow$Condition %in% levelList[[i]],]
  dataNow$Condition <- factor(dataNow$Condition, 
                              levels = levelList[[i]])
  if (max(dataNow$PHATE_1) < 0) {
    limsNow <- c((1.2*min(dataNow$PHATE_1)), 
      .8*max(dataNow$PHATE_1))
  } else {
    limsNow <- c((1.2*min(dataNow$PHATE_1)), 
                 1.3*max(dataNow$PHATE_1))
  }
  pltNow <- ggplot(data = dataNow, aes_string(x = "Condition",
                                              y = "PHATE_1")) +
    geom_boxplot() +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right", axis.text.x = element_text(size = 18)) +
    ylab("PHATE_1 Score") +
    scale_fill_manual(values=colors) +
    scale_y_continuous(expand = c(0,0), limits = limsNow) +
    labs(title = GSENow) + rremove('xlab') +
    theme(axis.text.y = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    rremove("legend") + rotate_x_text(45) +
    stat_compare_means(method = "t.test", label = "p.signif",
                       size = 7.5, hide.ns = FALSE,
                       method.args = list("alternative" = "greater"),
                       comparisons = compListList[[i]])
  pltList[[i]] <- pltNow
}
ga <- ggarrange(plotlist = pltList, ncol = 3, nrow = 2, align = "hv")
ggsave(ga, height = 16, width = 19, filename = "Figures_v2/PHATE_1_EWSFLI1sh.png")
ggsave(ga, height = 16, width = 19, filename = "Figures_v2/PHATE_1_EWSFLI1sh.pdf")

# Compare PHATE_1 for other interventions
GSEToPlot <- unique(phateEWS$SeriesGSE)
GSEToPlot <- c("GSE92739",
               "GSE98785", 
               "GSE118871",
               "GSE98786"
               )
dataNow <- phateEWS[phateEWS$SeriesGSE %in% GSEToPlot,]
pltList <- list()
compListList <- list(
  list(c("A673 siCTR", "A673 siMRTFA"),
       c("A673 siCTR", "A673 siMRTFB"),
       c("A673 siCTR", "A673 siMRTFA-siMRTFB"),
         c("A673 siCTR", "A673 siTEAD")),
  # list(c("A673 siCTR-shEF1", "A673 siMRTFA-shEF1"),
  #      c("A673 siCTR-shEF1", "A673 siMRTFB-shEF1"),
  #      c("A673 siCTR-shEF1", "A673 siMRTFA-siMRTFB-shEF1"),
  #      c("A673 siCTR-shEF1", "A673 siTEAD-shEF1")),
  rev(list(c("A673 DMSO", "A673 SP2509"),
       c("TC32 DMSO", "TC32 SP2509"),
       c("EWS502 DMSO", "EWS502 SP2509"),
       c("TC71 DMSO", "TC71 SP2509"))),
  list(c("A673 DMSO", "A673 SP2509"),
       c("A673 DMSO", "A673-SP2509R DMSO"),
       c("A673-SP2509R DMSO", "A673-SP2509R SP2509")),
  list(c("A673 iLuc", "A673 iLSD1")))
levelList <- list(
  c("A673 siCTR", "A673 siMRTFA", "A673 siMRTFB", "A673 siMRTFA-siMRTFB", "A673 siTEAD"),
  # c("A673 siCTR-shEF1", "A673 siMRTFA-shEF1", "A673 siMRTFB-shEF1", "A673 siMRTFA-siMRTFB-shEF1", "A673 siTEAD-shEF1"),
  NULL,
  c("A673 DMSO", "A673 SP2509", "A673-SP2509R DMSO", "A673-SP2509R SP2509"), NULL)
for (i in 1:length(GSEToPlot)) {
  GSENow <- GSEToPlot[i]
  dataNow <- phateEWS[phateEWS$SeriesGSE == GSENow,]
  dataNow <- dataNow[! is.na(dataNow$PHATE_1),]
  dataNow$Condition <- paste0(dataNow$Cell, " ", dataNow$Condition)
  if (i == 2) {
    dataNow$Condition <- gsub(dataNow$Condition, pattern = "HCI", replacement = "SP")
  }
  if (i == 3) {
    dataNow$Condition <- gsub(dataNow$Condition, pattern = "-0mo", replacement = "")
  }
  if(! is.null(levelList[[i]])) {
    levNow <- levelList[[i]]
  } else {
    levNow <- unique(dataNow$Condition)
  }
  dataNow <- dataNow[dataNow$Condition %in% levNow,]
  dataNow$Condition <- factor(dataNow$Condition, 
                              levels = levNow)
  if (max(dataNow$PHATE_1) < 0) {
    limsNow <- c((1.2*min(dataNow$PHATE_1)), 
                 .6*max(dataNow$PHATE_1))
  } else {
    limsNow <- c((1.2*min(dataNow$PHATE_1)), 
                 1.3*max(dataNow$PHATE_1))
  }
  if (i == 1) {
    limsNow <- c(-.0025, .001)
  } else if (i == 2) {
    limsNow <- c(-.006, 0.0045)
  } else if (i == 4) {
    limsNow <- c(-.004, .0008)
  } else if (i == 3) {
    limsNow <- c(-.004, .002)
  }
  pltNow <- ggplot(data = dataNow, aes_string(x = "Condition",
                                              y = "PHATE_1")) +
    geom_boxplot() +
    theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right", axis.text.x = element_text(size = 18)) +
    ylab("PHATE_1 Score") +
    scale_fill_manual(values=colors) +
    scale_y_continuous(expand = c(0,0), limits = limsNow) +
    labs(title = GSENow) + rremove('xlab') +
    theme(axis.text.y = element_text(size = 18),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    rremove("legend") + rotate_x_text(45) +
    stat_compare_means(method = "t.test", label = "p.signif",
                       size = 5, hide.ns = FALSE, 
                       comparisons = compListList[[i]])
  pltList[[i]] <- pltNow
}
ga <- ggarrange(plotlist = pltList, ncol = 2, nrow = 2, align = "hv")
ggsave(ga, height = 16, width = 14, filename = "Figures_v2/PHATE_1_Interventions.png")
ggsave(ga, height = 16, width = 14, filename = "Figures_v2/PHATE_1_Interventions.pdf")




# Howarth et al 2014
GSENow <- "GSE60949"
dataNow <- phateEWS[phateEWS$SeriesGSE == GSENow,]
dataNow <- dataNow[! is.na(dataNow$PHATE_1),]
pal <- brewer.pal(3, 'Greys')
dataNow$Condition <- factor(dataNow$Condition, levels = c("shCTR", "shEF1", "shEWSAT1"))
pltNow <- ggplot(data = dataNow, aes_string(x = "Condition",
                                            y = "PHATE_1",
                                            fill = "Condition")) +
  geom_boxplot() +
  scale_fill_manual(values = pal) +
  theme_pubr(border = T, base_size = 22, legend = "right") +
  ylab("PHATE_1 position") +
  labs(title = "Howarth et al. 2014") +
  rremove('xlab') +
  scale_y_continuous(limits = c(-0.005, 0.012)) +
  theme(axis.text.y = element_text(size = 18), 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  stat_compare_means(comparisons = list(c( "shEF1", "shCTR"),
                                        c("shEWSAT1", "shCTR")), 
                     bracket.size = .5, method.args = list("alternative" = "greater"),
                     method = "t.test", size = 8, label = "p.signif")
ggsave(pltNow, filename = "Figures_v2/Howarth2014.png", height = 7.5,
       width = 11)

# Pishas et al 2018
GSENow <- "GSE98785"
dataNow <- phateEWS[phateEWS$SeriesGSE == GSENow,]
dataNow <- dataNow[! is.na(dataNow$PHATE_1),]
pal <- brewer.pal(3, 'Greys')
dataNow$Condition[dataNow$Condition == "HCI2509"] <- "SP2509"
dataNow$Condition <- factor(dataNow$Condition,
                            levels = rev(c("SP2509", "DMSO")))
pltNow <- ggplot(data = dataNow, aes_string(x = "Cell",
                                            y = "PHATE_1",
                                            fill = "Condition")) +
  geom_boxplot() +
  scale_fill_manual(values = pal) +
  theme_pubr(border = T, base_size = 22, legend = "right") +
  ylab("PHATE_1 position") +
  labs(title = "Pishas et al. 2018") +
  rremove('xlab') +
  scale_y_continuous(limits = c(-0.006, 0.002)) +
  theme(axis.text.y = element_text(size = 18), 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(pltNow, filename = "Figures_v2/Pishas2018.png", height = 7.5,
       width = 11)


# Show tissue
pltDataNow <- pltDataNow[order(match(pltDataNow$TissueNow, rev(levels(pltDataNow$TissueNow)))),]
pltDataNow$Tissue <- pltDataNow$TissueNow
DimPlot2(data = pltDataNow, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                 color = "Tissue"), 
         plotName = "Fig1F", colorMap = colMapTissueNow)
DimPlot2(data = pltDataNow, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                 color = "Category"), colorMap = colMapCategory,
         plotName = "Fig1G")
DimPlot2(data = pltDataNonEWS, mapping = aes_string(x = "PHATE_1", y = "PHATE_2",
                                                    color = "Category"), colorMap = colMapCategory,
         plotName = "Fig1H")


# #loop through plots
# rename <- function(x){
#   if (x < 10) {
#     return(name <- paste('000',i,'plot.png',sep=''))
#   }
#   if (x < 100 && i >= 10) {
#     return(name <- paste('00',i,'plot.png', sep=''))
#   }
#   if (x >= 100) {
#     return(name <- paste('0', i,'plot.png', sep=''))
#   }
# }
# rename2 <- function(x){
#   if (x < 10) {
#     return(name <- paste('000',j,'plot.png',sep=''))
#   }
#   if (x < 100 && j >= 10) {
#     return(name <- paste('00',j,'plot.png', sep=''))
#   }
#   if (x >= 100) {
#     return(name <- paste('0', j,'plot.png', sep=''))
#   }
# }
# dir.create("animationFolder", showWarnings = FALSE)
# file.remove(list.files("animationFolder", full.names = TRUE))
# pltDataNow3D <- pltDataNow[pltDataNow$Category != "Mixed/other",]
# colMapCategory3D <- colMapCategory[c(-5)]
# frames <- 360
# frames2 <- 360
# for(i in 1:frames){
#   name <- paste0("animationFolder/", rename(i))
#   
#   #saves the plot as a .png file in the working directory
#   png(name, height = 5, width = 6, units = "in", res = 300)
#   par(mar = c(3,3,3,7))
#   scatter3D(x=pltDataNow3D$PHATE_1, scale = FALSE,
#             y=pltDataNow3D$PHATE_2, cex = .6,
#             z=pltDataNow3D$PHATE_3, pch = 16,
#             colvar = as.integer(pltDataNow3D$Category),
#             # colkey = list(rev(c("#C1C1C1", "#7CAE00", "#00BFC4", 
#             #                     "#C77CFF", "#F8766D")),
#             #               labels = c("setosa", "versicolor", "virginica")),
#             col = colMapCategory3D, 
#             colkey = list(length = .5, #side = 1, length = 1, 
#                           at = c(1.4, 2.2, 3, 3.8),
#                           labels = names(colMapCategory3D)),
#             # clab = "Category", 
#             phi = 30, d = 2,
#             theta = i, bty = "n")
#   dev.off()
#   if (i == frames) {
#     for(j in 1:frames){
#       name2 <- paste0("animationFolder/vert", rename2(j))
#       #saves the plot as a .png file in the working directory
#       png(name2, height = 5, width = 6, units = "in", res = 300)
#       par(mar = c(3,3,3,7))
#       scatter3D(x=pltDataNow3D$PHATE_1, scale = FALSE,
#                 y=pltDataNow3D$PHATE_2, cex = .6,
#                 z=pltDataNow3D$PHATE_3, pch = 16,
#                 colvar = as.integer(pltDataNow3D$Category),
#                 # colkey = list(rev(c("#C1C1C1", "#7CAE00", "#00BFC4", 
#                 #                     "#C77CFF", "#F8766D")),
#                 #               labels = c("setosa", "versicolor", "virginica")),
#                 col = colMapCategory3D, 
#                 colkey = list(length = .5, #side = 1, length = 1, 
#                               at = c(1.4, 2.2, 3, 3.8),
#                               labels = names(colMapCategory3D)),
#                 # clab = "Category", 
#                 phi = 30+j, d = 2,
#                 theta = i, bty = "n")
#       dev.off()
#     }
#   }
# }
# 
# my_command <- 'cd animationFolder && convert *.png -delay .1 -loop 0 ../PHATE_360.gif && cd ..'
# system(my_command)





cat("\n", timestamp2(), " Calculating gene correlations with PHATE and plotting...\n", sep = "")
if (! file.exists("Data/bulkRNASeq/trendyRes_NormalPHATE_WithEWSTrendy.rda")) {
  # DGEs and GSEA from PHATE analysis
  if (! "expression" %in% ls()) {
    load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  }
  pltDataNowEWS <- pltDataNow[pltDataNow$TissueFinal == "Ewing sarcoma",]
  countsNow <- expression[, colnames(expression) %in% 
                            pltDataNowEWS$samples]
  all(colnames(countsNow) == pltDataNowEWS$samples)
  Sizes <- MedianNorm(countsNow)
  countsNorm <- GetNormalizedMat(countsNow, Sizes)
  keep_rows <- rowSums(countsNorm > 0) > (length(colnames(countsNorm)) * .1)
  countsNorm <- countsNorm[keep_rows,]
  timeList <- list("PHATE_1" = pltDataNowEWS$PHATE_1, 
                   "PHATE_2" = pltDataNowEWS$PHATE_2)
  
  resListEWS <- lapply(timeList, FUN = function(timeNow) {
    doTrendyAnalysis(countsNorm, 
                     timeVec = timeNow, 
                     TERM2GENE = TERM2GENE)
  }) 
  save(resListEWS, file = "Data/bulkRNASeq/trendyRes_NormalPHATE_WithEWSTrendy.rda")
} else {
  load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithEWSTrendy.rda")
}


# Plot results with only EWS
pltListEWS <- lapply(names(resListEWS), 
                     FUN = doPhateCorrPlot, 
                     resSamples = resListEWS, intGenes = NULL)
g1 <- ggarrange(plotlist = pltListEWS, ncol = 3, nrow = 1)
ggsave(g1, height = 7.5, width = 24,
       filename = "Figures_v2/FigS1Mii.png", device = "png")


# Plot GSEA 
ranks <- resListEWS$PHATE_1$value
names(ranks) <- resListEWS$PHATE_1$geneName
GSEARes1 <- doAltGSEA(ranks, 
                      TERM2GENEList = TERM2GENEList)
results <- lapply(names(GSEARes1), FUN = function(nameNow) {
  resNow <- GSEARes1[[nameNow]][[2]]
  resNow$Group <- nameNow
  resNow
})
DF <- data.table::rbindlist(results)
eresSave <- DF[,c(2, 9, 3, 4, 5, 6, 7, 8)]
write.table(eresSave, file = "Tables/ewsGSEA_PHATE_1.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
plotList1 <- lapply(GSEARes1, FUN = function(x) {
  x$gseaPlot
})

ga1 <- ggarrange(plotlist = plotList1, align = "v",
                 nrow = 2, ncol = 2)
ggsave(ga1, filename = "Figures_v2/FigS1Miii.png",
       height = 12, width = 25)


if (! file.exists("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")) {
  # DGEs and GSEA from PHATE analysis
  if (! "expression" %in% ls()) {
    load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  }
  pltDataWithoutEWS <- pltDataNow[pltDataNow$TissueFinal != "Ewing sarcoma",]
  countsNow <- expression[, colnames(expression) %in% 
                            pltDataWithoutEWS$samples]
  all(colnames(countsNow) == pltDataWithoutEWS$samples)
  Sizes <- MedianNorm(countsNow)
  countsNorm <- GetNormalizedMat(countsNow, Sizes)
  keep_rows <- rowSums(countsNorm > 0) > (length(colnames(countsNorm)) * .1)
  countsNorm <- countsNorm[keep_rows,]
  timeList <- list("PHATE_1" = pltDataWithoutEWS$PHATE_1, 
                   "PHATE_2" = pltDataWithoutEWS$PHATE_2)
  
  resListWithoutEWS <- lapply(timeList, FUN = function(timeNow) {
    doTrendyAnalysis(countsNorm, 
                     timeVec = timeNow, 
                     TERM2GENE = TERM2GENE)
  }) 
  save(resListWithoutEWS, file = "Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
} else {
  load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
}

pltDF <- resListWithoutEWS$PHATE_1
pltDFEws <- resListEWS$PHATE_1
mergeFrame <- merge(pltDF[,c(-3)], pltDFEws[,c(-3)], by = "geneName", all = TRUE)
colnames(mergeFrame) <- c("geneName", "ewsCorrValue_PHATE_1", "nonEwsCorrValue_PHATE_1")
write.table(mergeFrame, file = "Tables/PHATE_1_correlationScores.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Plot results with only EWS
pltListWithoutEWS <- lapply(names(resListWithoutEWS), 
                            FUN = doPhateCorrPlot, 
                            resSamples = resListWithoutEWS,
                            intGenes = NULL)
g1 <- ggarrange(plotlist = pltListWithoutEWS, ncol = 3, nrow = 1)
ggsave(g1, height = 7.5, width = 24,
       filename = "Figures_v2/FigS1Miv.png", device = "png")

if (! "vsd" %in% ls()) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}
pal <- colorRampPalette(brewer.pal(9, "PuBu"))(255)
featurePlot2(pltData = pltDataNow, plotName = "FigS1Nii",
             x = "PHATE_1", y = "PHATE_2", pal = pal,
             countsNorm = vsd, features = "KLHL23")
pal2 <- colorRampPalette(brewer.pal(9, "OrRd"))(255)
featurePlot2(pltData = pltDataNow, plotName = "FigS1Niii",
             x = "PHATE_1", y = "PHATE_2",pal = pal2,
             countsNorm = vsd, features = "PARVA")
pal3 <- colorRampPalette(brewer.pal(9, "BuGn"))(255)
featurePlot2(pltData = pltDataNow, plotName = "FigS1Niv",
             x = "PHATE_1", y = "PHATE_2", pal = pal3,
             countsNorm = vsd, features = "MAP2")

# Plot GSEA 
ranks <- resListWithoutEWS$PHATE_1$value
names(ranks) <- resListWithoutEWS$PHATE_1$geneName
GSEARes1 <- doAltGSEA(ranks, 
                      TERM2GENEList = TERM2GENEList)
results <- lapply(names(GSEARes1), FUN = function(nameNow) {
  resNow <- GSEARes1[[nameNow]][[2]]
  resNow$Group <- nameNow
  resNow
})
DF <- data.table::rbindlist(results)
eresSaveNonEWS <- DF[,c(2, 9, 3, 4, 5, 6, 7, 8)]
write.table(eresSave, file = "Tables/nonEWSGSEA_PHATE_1.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

plotList1 <- lapply(GSEARes1, FUN = function(x) {
  x$gseaPlot
})

ga1 <- ggarrange(plotlist = plotList1, align = "v",
                 nrow = 2, ncol = 2)
ggsave(ga1, filename = "Figures_v2/FigS1Mv.png",
       height = 12, width = 25)

# Compare PHATE_1 for EWS and Non-EWS
pltDFEWS <- resListEWS$PHATE_1
pvec <- sapply(pltDFEWS$value, FUN = function(x, n) {
  # Convert sign r2 to r
  origVal <- x
  x <- sqrt(abs(x))*sign(origVal)
  dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
}, n = length(pltDFEWS$value))
pltDFEWS <- pltDFEWS[p.adjust(pvec, method = "BH") < .0001,]
colnames(pltDFEWS)[2] <- "EWS_value"
pltDFWithoutEWS <- resListWithoutEWS$PHATE_1
pvec <- sapply(pltDFWithoutEWS$value, FUN = function(x, n) {
  # Convert sign r2 to r
  origVal <- x
  x <- sqrt(abs(x))*sign(origVal)
  dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
}, n = length(pltDFWithoutEWS$value))
pltDFWithoutEWS <- pltDFWithoutEWS[p.adjust(pvec, method = "BH") < .0001,]
colnames(pltDFWithoutEWS)[2] <- "NONEWS_value"
mergePLT <- merge(pltDFEWS[,c(-3)], pltDFWithoutEWS[,c(-3)], by = "geneName")
lmNow <- lm(EWS_value ~ NONEWS_value, data = mergePLT)
summary(lmNow)
mergePLT$VAR <- abs(mergePLT$EWS_value - mergePLT$NONEWS_value)
mergePLT$switchGene <- FALSE
mergePLT$switchGene[sign(mergePLT$EWS_value) != sign(mergePLT$NONEWS_value)] <- TRUE
mergePLT$varGene <- FALSE
mergePLT$varGene[mergePLT$VAR > .35 & mergePLT$switchGene] <- TRUE
doMarksDF <- mergePLT[mergePLT$VAR > .5 & mergePLT$switchGene &
                        abs(mergePLT$EWS_value) > .15 & 
                        abs(mergePLT$NONEWS_value) > .15,]
g1 <- ggplot(mergePLT, mapping = aes_string(x = "EWS_value", y = "NONEWS_value",
                                            color = "varGene")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + 
  xlab("Ewing Sarcoma Phate_1 Score") +
  ylab("Developmental Context PHATE_1 Score") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,
  #                                                (max(mergePLT$NONEWS_value)*1.05))) +
  # scale_x_continuous(expand = c(0,0), limits = c(0, (max(mergePLT$)*1.05))) +
  scale_color_manual(values=c( "#C4BFBF", "#CC0909")) +
  labs(color = "Variant Gene") +
  geom_text_repel(data = doMarksDF,segment.colour = "black", 
                  color = "black", seed = 42, #nudge_x = .15,
                  # force = 2, point.padding = .5,
                  mapping = aes_string(x = "EWS_value", y = "NONEWS_value",
                                       label = "geneName"))
ggsave(g1, filename = "Figures_v2/FigS1Mvi.png",
       height = 7.5, width = 8)

## Compare GSEA results for PHATE_1 in EWS and Normal tissues
ranks <- resListEWS$PHATE_1$value
names(ranks) <- resListEWS$PHATE_1$geneName
GSEAResEWS <- doAltGSEA(ranks, 
                        TERM2GENEList = TERM2GENEList)
ranks <- resListWithoutEWS$PHATE_1$value
names(ranks) <- resListWithoutEWS$PHATE_1$geneName
GSEAResWithoutEWS <- doAltGSEA(ranks, 
                               TERM2GENEList = TERM2GENEList)

# Compare GO
GOEWS <- GSEAResEWS$CGP$eres
GONOEWS <- GSEAResWithoutEWS$CGP$eres
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES < 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES > 0])
ol <- calculate.overlap(oLap)
rmVec <- c(ol$a3)
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES > 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES < 0])
ol <- calculate.overlap(oLap)
rmVec <- unique(c(rmVec, ol$a3))
# Forcing non-overlap of opposing terms
GOEWS$ID[GOEWS$ID %in% rmVec] <- paste0(GOEWS$ID[GOEWS$ID %in% rmVec], "2")
GONOEWS$ID[GONOEWS$ID %in% rmVec] <- paste0(GONOEWS$ID[GONOEWS$ID %in% rmVec], "3")
# Now no terms going in opposite directions...
oLap <- list("Ewing sarcoma" = GOEWS$ID,
             "Developmental context" = GONOEWS$ID)
calculate.overlap.and.pvalue(oLap[[1]], oLap[[2]], lower.tail = FALSE,
                             total.size = length(unique(TERM2GENE_CGP$gs_name)))
gd <- venn.diagram(oLap, filename = NULL,cat.dist= c(.25, .25),
                   fill = c("Firebrick", "skyblue"))
plot.new()
dev.off()
grid.draw(gd)
gd <- grid.grab()
plt <- as.ggplot(gd)
dev.off()
ggsave(plt, filename = "Figures_v2/FigS1Mvii.png", height = 4, width = 6)
ggsave(plt, filename = "Figures_v2/FigS1Mvii.pdf", height = 4, width = 6)
# Shared Biological processes
GOEWS$Type <- "Ewing sarcoma"
topN <- 4
charLimit <- 50
eres <- GOEWS
eres <- eres[eres$ID %in% GONOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt1 <- rbind(eresUp, eresDn)
GONOEWS$Type <- "Developmental context"
eres <- GONOEWS
eres <- eres[eres$ID %in% GOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt2 <- rbind(eresUp, eresDn)
KeepTerms <- unique(c(eresPlt1$ID, eresPlt2$ID))
eres <- rbind(GOEWS, GONOEWS)
eres <- eres[eres$ID %in% KeepTerms,]
eres$Group <- "Under-expressed"
eres$Group[eres$NES > 0] <- "Over-expressed"
eres$NESEWING <- 0
eres$NESEWING[eres$Type == "Ewing sarcoma"] <- eres$NES[eres$Type == "Ewing sarcoma"]
eres <- eres[order(eres$NESEWING, decreasing = T),]
eres$ID <- fixStrings(eres$ID)
eres$ID <- gsub(eres$ID, pattern = "GO ", replacement = "")
eres$ID[nchar(eres$ID) > charLimit] <- paste0(substr(eres$ID[nchar(eres$ID) > charLimit], 1, (charLimit-3)), "...")
# eresPlt <- eresPlt[! duplicated(eresPlt$ID),]
eres$ID <- factor(eres$ID, levels = rev(unique(eres$ID[eres$NESEWING != 0])))
g1 <- ggplot(data = eres, 
             mapping = aes_string(x = "ID", 
                                  y = "NES", 
                                  fill = "Type")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Normalized Enrichment Score") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  scale_fill_manual(values=c("skyblue","firebrick")) +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(expand = c(0,0), limits = c((min(eres$NES)*1.1),(max(eres$NES)*1.1))) 
ggsave(g1, filename = "Figures_v2/FigS1Mviii.png",
       height = 7, width = 12)


# Compare EWS
GOEWS <- GSEAResEWS$`Ewing sarcoma`$eres
GONOEWS <- GSEAResWithoutEWS$`Ewing sarcoma`$eres
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES < 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES > 0])
ol <- calculate.overlap(oLap)
print(ol$a3)
rmVec <- c(ol$a3)
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES > 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES < 0])
ol <- calculate.overlap(oLap)
print(ol$a3)
rmVec <- unique(c(rmVec, ol$a3))
# Only 1 term was opposite -- removing
GOEWS$ID[GOEWS$ID %in% rmVec] <- paste0(GOEWS$ID[GOEWS$ID %in% rmVec], "2")
GONOEWS$ID[GONOEWS$ID %in% rmVec] <- paste0(GONOEWS$ID[GONOEWS$ID %in% rmVec], "3")
# Now no terms going in opposite directions...
oLap <- list("Ewing sarcoma" = GOEWS$ID,
             "Developmental context" = GONOEWS$ID)
gd <- venn.diagram(oLap, filename = NULL,cat.dist= c(.25, .25),
                   fill = c("Firebrick", "skyblue"))
calculate.overlap.and.pvalue(oLap[[1]], oLap[[2]], lower.tail = F,
                             total.size = length(unique(TERM2GENE_EWS$gs_name)))
plot.new()
dev.off()
grid.draw(gd)
gd <- grid.grab()
plt <- as.ggplot(gd)
dev.off()
ggsave(plt, filename = "Figures_v2/FigS1Mvix.png", height = 4, width = 6)
ggsave(plt, filename = "Figures_v2/FigS1Mvix.pdf", height = 4, width = 6)
# Shared Biological processes
GOEWS$Type <- "Ewing sarcoma"
topN <- 4
charLimit <- 50
eres <- GOEWS
eres <- eres[eres$ID %in% GONOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresUp <- eresUp[eresUp$NES > 0,]
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresDn <- eresDn[eresDn$NES < 0,]
eresPlt1 <- rbind(eresUp, eresDn)
GONOEWS$Type <- "Developmental context"
eres <- GONOEWS
eres <- eres[eres$ID %in% GOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresUp <- eresUp[eresUp$NES > 0,]
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresDn <- eresDn[eresDn$NES < 0,]
eresPlt2 <- rbind(eresUp, eresDn)
KeepTerms <- unique(c(eresPlt1$ID, eresPlt2$ID))
eres <- rbind(GOEWS, GONOEWS)
eres <- eres[eres$ID %in% KeepTerms,]
eres$Group <- "Under-expressed"
eres$Group[eres$NES > 0] <- "Over-expressed"
eres$NESEWING <- 0
eres$NESEWING[eres$Type == "Ewing sarcoma"] <- eres$NES[eres$Type == "Ewing sarcoma"]
eres <- eres[order(eres$NESEWING, decreasing = T),]
eres$ID <- fixStrings(eres$ID)
eres$ID <- gsub(eres$ID, pattern = "GO ", replacement = "")
eres$ID[nchar(eres$ID) > charLimit] <- paste0(substr(eres$ID[nchar(eres$ID) > charLimit], 1, (charLimit-3)), "...")
# eresPlt <- eresPlt[! duplicated(eresPlt$ID),]
eres$ID <- factor(eres$ID, levels = rev(unique(eres$ID[eres$NESEWING != 0])))
g1 <- ggplot(data = eres, 
             mapping = aes_string(x = "ID", 
                                  y = "NES", 
                                  fill = "Type")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Normalized Enrichment Score") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  scale_fill_manual(values=c("skyblue","firebrick")) +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(expand = c(0,0), limits = c((min(eres$NES)*1.1),(max(eres$NES)*1.1))) 

ggsave(g1, filename = "Figures_v2/FigS1Mx.png",
       height = 7, width = 12)


#####################################################################################
########################## Hallmarks of EWS-FLI1 Permissibility #####################
#####################################################################################

# Step 1: Transcriptional phenotypes associated with EWS-FLI1 expression
ef1UpGenes <- TERM2GENE_EWS$gene_symbol[grep(x = TERM2GENE_EWS$gs_name, 
                                             pattern = "_UP|ACTIVATED|HIGH")]
ef1UpPaths <- enricher(unique(ef1UpGenes), universe = unique(TERM2GENE_EWS$gene_symbol),
                       TERM2GENE = TERM2GENE)
ef1UpPaths <- as.data.frame(ef1UpPaths)
ef1UpPaths$Group <- "EWS-FLI1 transcriptome"

# Step 2: Transcriptional phenotypes associated with PHATE_1-low vs PHATE_1-high
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
rankDF <- resListWithoutEWS$PHATE_1
PHATE1DnGenes <- rankDF$geneName[rankDF$value < -.2]
PHATE1DnPaths <- enricher(PHATE1DnGenes,
                          TERM2GENE = TERM2GENE,
                          universe = unique(rankDF$geneName))
PHATE1DnPaths <- as.data.frame(PHATE1DnPaths)

# # Step 3: Transcriptional phenotypes associated with PHATE_1-high vs PHATE_1-low
# PHATE1UpGenes <- rankDF$geneName[rankDF$value > .2]
# PHATE1UpPaths <- enricher(PHATE1UpGenes, TERM2GENE = TERM2GENE)
# PHATE1UpPaths <- as.data.frame(PHATE1UpPaths)

# Step 4: Compare PHATE_1-low and EWS-FLI1 pathways
oList <- list("EWS-FLI1 activated" = ef1UpPaths$ID,
              "PHATE_1-low" = PHATE1DnPaths$ID)
calculate.overlap.and.pvalue(oList[[1]], oList[[2]], lower.tail = F,
                             total.size = length(unique(TERM2GENE$gs_name)))
gd <- venn.diagram(oList, filename = NULL,cat.dist= c(.25, .25),
                   fill = c("#F8766D", "#C77CFF"))
plot.new()
dev.off()
grid.draw(gd)
gd <- grid.grab()
plt <- as.ggplot(gd)
dev.off()
ggsave(plt, filename = "Figures_v2/Venn_EWSFL1_Activated_PHATE1_Low.pdf", height = 4, width = 6)
sharedPathsDn <- calculate.overlap(oList)
sharedPathsDn <- sharedPathsDn$a3

# # Step 5: Compare PHATE_1-high and EWS-FLI1 pathways
# oList <- list("EWS-FLI1 activated" = ef1UpPaths$ID,
#               "PHATE_1-high" = PHATE1UpPaths$ID)
# calculate.overlap.and.pvalue(oList[[1]], oList[[2]], lower.tail = F,
#                              total.size = length(unique(TERM2GENE$gs_name)))
# gd <- venn.diagram(oList, filename = NULL,cat.dist= c(.25, .25),
#                    fill = c("#F8766D", "#C77CFF"))
# plot.new()
# dev.off()
# grid.draw(gd)
# gd <- grid.grab()
# plt <- as.ggplot(gd)
# dev.off()
# ggsave(plt, filename = "Figures_v2/Venn_EWSFL1_Activated_PHATE1_High.pdf", height = 4, width = 6)
# sharedPathsUp <- calculate.overlap(oList)
# sharedPathsUp <- sharedPathsUp$a3

# Step 6: combined pathway results
# sharedPaths <- unique(c(sharedPathsUp, sharedPathsDn))

# Bar chart of shared pathways
sharedPathsDF <- data.frame(
  "Shared_Path" = sharedPathsDn,
  "Category" = rep("Other", length(sharedPathsDn)), 
  stringsAsFactors = FALSE
)
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
     pattern = "E2F|Cycle|G2M|PROLIF|_PHASE|MITOT|G2|KINETOCH|CENTRO|REPLICAT|TELOM|CHROMOSOME|CHROMATID|DIVID|CHROMOSOMAL")
     ] <- "Cell cycle"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "_HIV_")
                       ] <- "Other"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "CHROMATIN|EZH2|HISTONE|METHYL|SWI")
                       ] <- "Chromatin remodeling"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "ERROR|REPAIR|ANNEAL|HOMOLOGOUS_DNA_PAIRING|_ATR|INTEGRITY|DAMAGE|RESOLUTION|ATM|BRCA|RECOMB|FANCO")
                       ] <- "DNA repair"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "CANCER|OMA_|OMA$|METAST")
                       ] <- "Cancer"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "ES_1|PLURI|STEM")
                       ] <- "Stemness"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "HEAT")
                       ] <- "Heat shock"
sharedPathsDF$Category[grep(sharedPathsDF$Shared_Path, ignore.case = TRUE,
                            pattern = "_RNA|SPLIC|_MRNA")
                       ] <- "RNA processing"
numberPaths <- as.data.frame(table(sharedPathsDF$Category))
numberPaths <- numberPaths[! numberPaths$Var1 == "Other",]
colnames(numberPaths) <- c("Category", "Number")

numberPaths <- numberPaths[order(numberPaths$Number, decreasing = T),]
numberPaths$Category <- factor(numberPaths$Category, levels = rev(numberPaths$Category))
g1 <- ggplot(data = numberPaths, 
             mapping = aes_string(x = "Category", 
                                  y = "Number")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Number shared gene sets") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures_v2/EWSFLI1_PHATE1Dn_Shared_PathEnrich.pdf",
       height = 7.5, width = 8)

# DNA Repair pathways
g1 <- ggplot(data = numberPaths, 
             mapping = aes_string(x = "Category", 
                                  y = "Number")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Number shared gene sets") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures_v2/EWSFLI1_PHATE1Dn_Shared_PathEnrich.pdf", height = 7.5, width = 8)
load("Data/bulkRNASeq/fullUMAPData.rda")
# Classify tissue types
tissueDictNow <- jsonlite::read_json("Data/bulkRNASeq/tissueDictionary.json")
pltDataToClass <- unique(pltData[,c(1:3)])
pltDataToClass <- categorizeMetaData(pltDataToClass,
                                     cols = c("tissue", "TissueType"), 
                                     dictionary = tissueDictNow[c(1:6,8:15,17:31, 33:39)])
rSums <- rowSums(pltDataToClass[,c(4:39)])
toClass <- pltDataToClass[rSums == 1,]
classes <- c()
toClass$Tissue <- "None"
possibles <- colnames(toClass[,c(4:39)])
resVec <- c()
for (i in 1:length(toClass$samples)) {
  resVec <- c(resVec, possibles[which(toClass[i,c(4:39)] == 1)])
}
toClass$Tissue <- resVec
table(toClass$Tissue)
classMap <- toClass[,c(1, 2, 3, 40)]
pltData <- pltData[,c(-13)]
catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "stem-like"] <- "Other stem-cells"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "MSCs"] <- "MSCs"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "HSCs"] <- "HSCs"
table(catClass$Tissue)
pltData$TissueFinal <- catClass$Tissue
pltData <- pltData[pltData$TissueFinal != "Non-Ewing tumor",]
pltData$TissueFinal[is.na(pltData$TissueFinal)] <- "Mixed/other"
if (! 'vsd' %in% ls()) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}
vsd <- vsd[, colnames(vsd) %in% pltData$samples]
pltData <- pltData[pltData$samples %in% colnames(vsd),]
all(pltData$samples == colnames(vsd))
pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                 "iPSCs", "hESCs")] <- 
  str_to_sentence(pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                                   "iPSCs", "hESCs")])

all(pltData$samples == colnames(vsd))
pltData$DDR <- colMeans(vsd[c("FANCA", "FANCD2", 
                      "FANCI", "FEN1"),])
pltData$Tissue <- pltData$TissueFinal
# Save expression info for BLM FEN1 etc
colnames(pltData)
# expSave <- pltData[, c(1, 2, 19, 14:17)]
# colnames(expSave)[c(1:2)] <- c("Samples", "Series")
# write.table(expSave, file = "Tables/RepStressGeneExp.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
meds <- sapply(unique(pltData$Tissue), FUN = function(Tissue) {
  median(pltData$DDR[pltData$Tissue == Tissue])
})
names(meds) <- unique(pltData$Tissue)
meds <- meds[order(meds)]
meds <- c(meds[! names(meds) == "Ewing sarcoma"], meds[ names(meds) == "Ewing sarcoma"])
pltData <- pltData[order(match(pltData$Tissue, names(meds))),]
pltData$Tissue <- factor(pltData$Tissue, levels = unique(pltData$Tissue))
g1 <- ggplot(pltData, mapping = aes_string(y = "DDR",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "DNA Repair Genes") +
  rremove("xlab") + rotate_x_text(angle = 45)
g1

ggsave(g1, filename = "Figures_v2/FigS1X_repStressGenes_1.pdf",
       height = 7.5, width = 12)
g2 <- ggplot(pltData, mapping = aes_string(y = "FANCI",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCI") +
  rremove("xlab") + rotate_x_text(angle = 45)
g2
ggsave(g2, filename = "Figures_v2/FigS1X_repStressGenes_2.pdf",
       height = 7.5, width = 12)
g3 <- ggplot(pltData, mapping = aes_string(y = "FANCA",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCA") +
  rremove("xlab") + rotate_x_text(angle = 45)
g3
ggsave(g3, filename = "Figures_v2/FigS1X_repStressGenes_3.pdf",
       height = 7.5, width = 12)
g4 <- ggplot(pltData, mapping = aes_string(y = "FANCD2",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCD2") +
  rremove("xlab") + rotate_x_text(angle = 45)
g4
ggsave(g4, filename = "Figures_v2/FigS1X_repStressGenes_4.pdf",
       height = 7.5, width = 12)

ga <- ggarrange(g1, g2, g3, g4, ncol = 1, nrow = 4)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.png",
       height = 30, width = 12)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.pdf",
       height = 30, width = 12)





# Replication stress pathways
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithEWSTrendy.rda")
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
TERM2GENE_RS <- TERM2GENE[grep(TERM2GENE$gs_name, ignore.case = TRUE,
                               pattern = "replication_stress|dna_damage|repair|fork"),]
rankDF <- resListEWS$PHATE_1
ranks <- rankDF$value
names(ranks) <- rankDF$geneName
resRS <- doAltGSEA(ranks, nCharLimit = 60,
                   TERM2GENEList = list("Replication stress gene sets" = TERM2GENE_RS))
TERM2GENEList[["Replication stress gene sets"]] <- TERM2GENE_RS
ranks <- resListEWS$PHATE_1$value
names(ranks) <- resListEWS$PHATE_1$geneName
resNonRS <- doAltGSEA(ranks, 
                      TERM2GENEList = TERM2GENEList)
results <- lapply(names(resNonRS), FUN = function(nameNow) {
  resNow <- resNonRS[[nameNow]][[2]]
  resNow$Group <- nameNow
  resNow
})
DF <- data.table::rbindlist(results)
eresSaveEWS <- DF[,c(2, 9, 3, 4, 5, 6, 7, 8)]
write.table(eresSaveEWS, file = "Tables/ewsGSEA_PHATE_1.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
ranks <- resListWithoutEWS$PHATE_1$value
names(ranks) <- resListWithoutEWS$PHATE_1$geneName
resNonRS <- doAltGSEA(ranks, 
                      TERM2GENEList = TERM2GENEList)
results <- lapply(names(resNonRS), FUN = function(nameNow) {
  resNow <- resNonRS[[nameNow]][[2]]
  resNow$Group <- nameNow
  resNow
})
DF <- data.table::rbindlist(results)
eresSaveNonEWS <- DF[,c(2, 9, 3, 4, 5, 6, 7, 8)]
write.table(eresSaveNonEWS, file = "Tables/nonEwsGSEA_PHATE_1.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

plt <- resRS$`Replication stress gene sets`$gseaPlot
ggsave(plt, filename = "Figures_v2/ReplicationStress_EWS_GSEA.png", height = 6, width = 12)
rankDF <- resListWithoutEWS$PHATE_1
ranks <- rankDF$value
names(ranks) <- rankDF$geneName
resRS2 <- doAltGSEA(ranks, nCharLimit = 60,
                    TERM2GENEList = list("Replication stress gene sets" = TERM2GENE_RS))
plt <- resRS2$`Replication stress gene sets`$gseaPlot
ggsave(plt, filename = "Figures_v2/ReplicationStress_NOEWS_GSEA.png", height = 6, width = 12)

# Compare
GOEWS <- resRS$`Replication stress gene sets`$eres
GONOEWS <- resRS2$`Replication stress gene sets`$eres
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES < 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES > 0])
ol <- calculate.overlap(oLap)
rmVec <- c(ol$a3)
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES > 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES < 0])
ol <- calculate.overlap(oLap)
rmVec <- unique(c(rmVec, ol$a3))
# No terms going in opposite directions...
oLap <- list("Ewing sarcoma" = GOEWS$ID,
             "Developmental context" = GONOEWS$ID)
ol <- calculate.overlap(oLap)
# Calculate overlap p values
calculate.overlap.and.pvalue(oLap[[1]], oLap[[2]],
                             length(unique(TERM2GENE_CGP$gs_name)),
                             lower.tail = FALSE)


gd <- venn.diagram(oLap, filename = NULL,cat.dist= c(.25, .25),
                   fill = c("Firebrick", "skyblue"))
plot.new()
dev.off()
grid.draw(gd)
gd <- grid.grab()
plt <- as.ggplot(gd)
dev.off()
ggsave(plt, filename = "Figures_v2/Fig4A.png", height = 4, width = 6)
# Shared Biological processes
GOEWS$Type <- "Ewing sarcoma"
topN <- 5
charLimit <- 75
eres <- GOEWS
eres <- eres[eres$ID %in% GONOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt1 <- rbind(eresUp, eresDn)
GONOEWS$Type <- "Developmental context"
eres <- GONOEWS
eres <- eres[eres$ID %in% GOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt2 <- rbind(eresUp, eresDn)
KeepTerms <- unique(c(eresPlt1$ID, eresPlt2$ID))
eres <- rbind(GOEWS, GONOEWS)
eres <- eres[eres$ID %in% KeepTerms,]
eres$Group <- "Under-expressed"
eres$Group[eres$NES > 0] <- "Over-expressed"
eres$NESEWING <- 0
eres$NESEWING[eres$Type == "Ewing sarcoma"] <- eres$NES[eres$Type == "Ewing sarcoma"]
eres <- eres[order(eres$NESEWING, decreasing = T),]
eres$ID <- fixStrings(eres$ID)
eres$ID <- gsub(eres$ID, pattern = "GO ", replacement = "")
eres$ID[nchar(eres$ID) > charLimit] <- paste0(substr(eres$ID[nchar(eres$ID) > charLimit], 1, (charLimit-3)), "...")
# eresPlt <- eresPlt[! duplicated(eresPlt$ID),]
eres$ID <- factor(eres$ID, levels = rev(unique(eres$ID[eres$NESEWING != 0])))
g1 <- ggplot(data = eres, 
             mapping = aes_string(x = "ID", 
                                  y = "NES", 
                                  fill = "Type")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Normalized Enrichment Score") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  scale_fill_manual(values=c("skyblue","firebrick")) +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(expand = c(0,0), limits = c((min(eres$NES)*1.1),(max(eres$NES)*1.1))) 

ggsave(g1, filename = "Figures_v2/Fig4B.png",
       height = 7, width = 12)

load("Data/bulkRNASeq/fullUMAPData.rda")
# Classify tissue types
tissueDictNow <- jsonlite::read_json("Data/bulkRNASeq/tissueDictionary.json")
pltDataToClass <- unique(pltData[,c(1:3)])
pltDataToClass <- categorizeMetaData(pltDataToClass,
                                     cols = c("tissue", "TissueType"), 
                                     dictionary = tissueDictNow[c(1:6,8:15,17:31, 33:39)])
rSums <- rowSums(pltDataToClass[,c(4:39)])
toClass <- pltDataToClass[rSums == 1,]
classes <- c()
toClass$Tissue <- "None"
possibles <- colnames(toClass[,c(4:39)])
resVec <- c()
for (i in 1:length(toClass$samples)) {
  resVec <- c(resVec, possibles[which(toClass[i,c(4:39)] == 1)])
}
toClass$Tissue <- resVec
table(toClass$Tissue)
classMap <- toClass[,c(1, 2, 3, 40)]
pltData <- pltData[,c(-13)]
catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "stem-like"] <- "Other stem-cells"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "MSCs"] <- "MSCs"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "HSCs"] <- "HSCs"
table(catClass$Tissue)
pltData$TissueFinal <- catClass$Tissue
pltData <- pltData[pltData$TissueFinal != "Non-Ewing tumor",]
pltData$TissueFinal[is.na(pltData$TissueFinal)] <- "Mixed/other"
if (! 'vsd' %in% ls()) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}
vsd <- vsd[, colnames(vsd) %in% pltData$samples]
pltData <- pltData[pltData$samples %in% colnames(vsd),]
all(pltData$samples == colnames(vsd))
pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                 "iPSCs", "hESCs")] <- 
  str_to_sentence(pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                                   "iPSCs", "hESCs")])

all(pltData$samples == colnames(vsd))
pltData$BLM <- vsd["BLM",]
pltData$FEN1 <- vsd["FEN1",]
pltData$FANCA <- vsd["FANCA",]
pltData$FANCI <- vsd["FANCI",]
pltData$FANCD2 <- vsd["FANCD2",]
pltData$BRCA1 <- vsd["BRCA1",]
pltData$Tissue <- pltData$TissueFinal

# Save expression info for BLM FEN1 etc
colnames(pltData)
expSave <- pltData[, c(1, 2, 19, 14:17)]
colnames(expSave)[c(1:2)] <- c("Samples", "Series")
write.table(expSave, file = "Tables/RepStressGeneExp.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


meds <- sapply(unique(pltData$Tissue), FUN = function(Tissue) {
  mean(c(median(pltData$FEN1[pltData$Tissue == Tissue]),
         median(pltData$FANCA[pltData$Tissue == Tissue]),
         median(pltData$FANCI[pltData$Tissue == Tissue]),
         median(pltData$FANCD2[pltData$Tissue == Tissue])))
})
names(meds) <- unique(pltData$Tissue)
meds <- meds[order(meds)]
meds <- c(meds[! names(meds) == "Ewing sarcoma"], meds[ names(meds) == "Ewing sarcoma"])
pltData <- pltData[order(match(pltData$Tissue, names(meds))),]
pltData$Tissue <- factor(pltData$Tissue, levels = unique(pltData$Tissue))
g1 <- ggplot(pltData, mapping = aes_string(y = "FEN1",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FEN1") +
  rremove("xlab") + rotate_x_text(angle = 45)
g1
ggsave(g1, filename = "Figures_v2/FigS1X_repStressGenes_1.pdf",
       height = 7.5, width = 12)
g2 <- ggplot(pltData, mapping = aes_string(y = "FANCI",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCI") +
  rremove("xlab") + rotate_x_text(angle = 45)
g2
ggsave(g2, filename = "Figures_v2/FigS1X_repStressGenes_2.pdf",
       height = 7.5, width = 12)
g3 <- ggplot(pltData, mapping = aes_string(y = "FANCA",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCA") +
  rremove("xlab") + rotate_x_text(angle = 45)
g3
ggsave(g3, filename = "Figures_v2/FigS1X_repStressGenes_3.pdf",
       height = 7.5, width = 12)
g4 <- ggplot(pltData, mapping = aes_string(y = "FANCD2",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCD2") +
  rremove("xlab") + rotate_x_text(angle = 45)
g4
ggsave(g4, filename = "Figures_v2/FigS1X_repStressGenes_4.pdf",
       height = 7.5, width = 12)

ga <- ggarrange(g1, g2, g3, g4, ncol = 1, nrow = 4)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.png",
       height = 30, width = 12)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.pdf",
       height = 30, width = 12)


# RNA processing pathways
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithEWSTrendy.rda")
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
TERM2GENE_RP <- TERM2GENE[grep(TERM2GENE$gs_name, ignore.case = TRUE,
                               pattern = "_RNA_|^RNA_|_RNA$|Splicing"),]
# TERM2GENE_RP <- TERM2GENE
TERM2GENEList[["RNA Processing"]] <- TERM2GENE_RP
rankDF <- resListEWS$PHATE_1
ranks <- rankDF$value
names(ranks) <- rankDF$geneName
resRS <- doAltGSEA(ranks, nCharLimit = 60,
                   TERM2GENEList = list("RNA Processing" = TERM2GENE_RP))
eresNow <- resRS$`RNA Processing`$eres
plt <- resRS$`RNA Processing`$gseaPlot
ggsave(plt, filename = "Figures_v2/RNAProcess_EWS_GSEA.png", height = 6, width = 12)
rankDF <- resListWithoutEWS$PHATE_1
ranks <- rankDF$value
names(ranks) <- rankDF$geneName
resRS2 <- doAltGSEA(ranks, nCharLimit = 60,
                    TERM2GENEList = list("RNA Processing" = TERM2GENE_RP))
plt <- resRS2$`RNA Processing`$gseaPlot
ggsave(plt, filename = "Figures_v2/RNAProcess_NOEWS_GSEA.png", height = 6, width = 12)

# Compare
GOEWS <- resRS$`RNA Processing`$eres
GONOEWS <- resRS2$`RNA Processing`$eres
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES < 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES > 0])
ol <- calculate.overlap(oLap)
rmVec <- c(ol$a3)
oLap <- list("Ewing sarcoma" = GOEWS$ID[GOEWS$NES > 0],
             "Developmental context" = GONOEWS$ID[GONOEWS$NES < 0])
ol <- calculate.overlap(oLap)
rmVec <- unique(c(rmVec, ol$a3))
# No terms going in opposite directions...
oLap <- list("Ewing sarcoma" = GOEWS$ID,
             "Developmental context" = GONOEWS$ID)
ol <- calculate.overlap(oLap)
# Calculate overlap p values
calculate.overlap.and.pvalue(oLap[[1]], oLap[[2]],
                             length(unique(TERM2GENE_CGP$gs_name)),
                             lower.tail = FALSE)


gd <- venn.diagram(oLap, filename = NULL,cat.dist= c(.25, .25),
                   fill = c("Firebrick", "skyblue"))
plot.new()
dev.off()
grid.draw(gd)
gd <- grid.grab()
plt <- as.ggplot(gd)
dev.off()
ggsave(plt, filename = "Figures_v2/RNAProcessing_Venn_Diagram.png", height = 4, width = 6)
# Shared Biological processes
GOEWS$Type <- "Ewing sarcoma"
topN <- 5
charLimit <- 75
eres <- GOEWS
eres <- eres[eres$ID %in% GONOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt1 <- rbind(eresUp, eresDn)
GONOEWS$Type <- "Developmental context"
eres <- GONOEWS
eres <- eres[eres$ID %in% GOEWS$ID,]
eresUp <- eres %>% top_n(n = topN, wt = NES)
eresDn <- eres %>% top_n(n = topN, wt = -NES)
eresPlt2 <- rbind(eresUp, eresDn)
KeepTerms <- unique(c(eresPlt1$ID, eresPlt2$ID))
eres <- rbind(GOEWS, GONOEWS)
eres <- eres[eres$ID %in% KeepTerms,]
eres$Group <- "Under-expressed"
eres$Group[eres$NES > 0] <- "Over-expressed"
eres$NESEWING <- 0
eres$NESEWING[eres$Type == "Ewing sarcoma"] <- eres$NES[eres$Type == "Ewing sarcoma"]
eres <- eres[order(eres$NESEWING, decreasing = T),]
eres$ID <- fixStrings(eres$ID)
eres$ID <- gsub(eres$ID, pattern = "GO ", replacement = "")
eres$ID[nchar(eres$ID) > charLimit] <- paste0(substr(eres$ID[nchar(eres$ID) > charLimit], 1, (charLimit-3)), "...")
# eresPlt <- eresPlt[! duplicated(eresPlt$ID),]
eres$ID <- factor(eres$ID, levels = rev(unique(eres$ID[eres$NESEWING != 0])))
g1 <- ggplot(data = eres, 
             mapping = aes_string(x = "ID", 
                                  y = "NES", 
                                  fill = "Type")) +
  geom_bar(stat = "identity", position = position_dodge(1)) + 
  ylab("Normalized Enrichment Score") +
  ggpubr::rotate() +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.position="top") + rremove("ylab") +
  scale_fill_manual(values=c("skyblue","firebrick")) +
  labs(fill = NULL) +
  theme(axis.text.y = element_text(size = 18),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_continuous(expand = c(0,0), limits = c((min(eres$NES)*1.1),(max(eres$NES)*1.1))) 
ggsave(g1, filename = "Figures_v2/RNA_Processing_Shared_GSEA.png",
       height = 7, width = 12)

load("Data/bulkRNASeq/fullUMAPData.rda")
# Classify tissue types
tissueDictNow <- jsonlite::read_json("Data/bulkRNASeq/tissueDictionary.json")
pltDataToClass <- unique(pltData[,c(1:3)])
pltDataToClass <- categorizeMetaData(pltDataToClass,
                                     cols = c("tissue", "TissueType"), 
                                     dictionary = tissueDictNow[c(1:6,8:15,17:31, 33:39)])
rSums <- rowSums(pltDataToClass[,c(4:39)])
toClass <- pltDataToClass[rSums == 1,]
classes <- c()
toClass$Tissue <- "None"
possibles <- colnames(toClass[,c(4:39)])
resVec <- c()
for (i in 1:length(toClass$samples)) {
  resVec <- c(resVec, possibles[which(toClass[i,c(4:39)] == 1)])
}
toClass$Tissue <- resVec
table(toClass$Tissue)
classMap <- toClass[,c(1, 2, 3, 40)]
pltData <- pltData[,c(-13)]
catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "stem-like"] <- "Other stem-cells"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "MSCs"] <- "MSCs"
catClass$Tissue[is.na(catClass$Tissue) & 
                  catClass$TissueType == "HSCs"] <- "HSCs"
table(catClass$Tissue)
pltData$TissueFinal <- catClass$Tissue
pltData <- pltData[pltData$TissueFinal != "Non-Ewing tumor",]
pltData$TissueFinal[is.na(pltData$TissueFinal)] <- "Mixed/other"
if (! 'vsd' %in% ls()) {
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}
vsd <- vsd[, colnames(vsd) %in% pltData$samples]
pltData <- pltData[pltData$samples %in% colnames(vsd),]
all(pltData$samples == colnames(vsd))
pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                 "iPSCs", "hESCs")] <- 
  str_to_sentence(pltData$TissueFinal[! pltData$TissueFinal %in% c("MSCs", "HSCs",
                                                                   "iPSCs", "hESCs")])

all(pltData$samples == colnames(vsd))
pltData$BLM <- vsd["BLM",]
pltData$FEN1 <- vsd["FEN1",]
pltData$FANCA <- vsd["FANCA",]
pltData$FANCI <- vsd["FANCI",]
pltData$FANCD2 <- vsd["FANCD2",]
pltData$BRCA1 <- vsd["BRCA1",]
pltData$Tissue <- pltData$TissueFinal

# Save expression info for BLM FEN1 etc
colnames(pltData)
expSave <- pltData[, c(1, 2, 19, 14:17)]
colnames(expSave)[c(1:2)] <- c("Samples", "Series")
write.table(expSave, file = "Tables/RepStressGeneExp.tsv", quote = FALSE, row.names = FALSE, sep = "\t")


meds <- sapply(unique(pltData$Tissue), FUN = function(Tissue) {
  mean(c(median(pltData$FEN1[pltData$Tissue == Tissue]),
         median(pltData$FANCA[pltData$Tissue == Tissue]),
         median(pltData$FANCI[pltData$Tissue == Tissue]),
         median(pltData$FANCD2[pltData$Tissue == Tissue])))
})
names(meds) <- unique(pltData$Tissue)
meds <- meds[order(meds)]
meds <- c(meds[! names(meds) == "Ewing sarcoma"], meds[ names(meds) == "Ewing sarcoma"])
pltData <- pltData[order(match(pltData$Tissue, names(meds))),]
pltData$Tissue <- factor(pltData$Tissue, levels = unique(pltData$Tissue))
g1 <- ggplot(pltData, mapping = aes_string(y = "FEN1",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FEN1") +
  rremove("xlab") + rotate_x_text(angle = 45)
g1
ggsave(g1, filename = "Figures_v2/FigS1X_repStressGenes_1.pdf",
       height = 7.5, width = 12)
g2 <- ggplot(pltData, mapping = aes_string(y = "FANCI",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCI") +
  rremove("xlab") + rotate_x_text(angle = 45)
g2
ggsave(g2, filename = "Figures_v2/FigS1X_repStressGenes_2.pdf",
       height = 7.5, width = 12)
g3 <- ggplot(pltData, mapping = aes_string(y = "FANCA",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCA") +
  rremove("xlab") + rotate_x_text(angle = 45)
g3
ggsave(g3, filename = "Figures_v2/FigS1X_repStressGenes_3.pdf",
       height = 7.5, width = 12)
g4 <- ggplot(pltData, mapping = aes_string(y = "FANCD2",
                                           x = "Tissue")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  ylab("Normalized read counts")  +
  labs(title = "FANCD2") +
  rremove("xlab") + rotate_x_text(angle = 45)
g4
ggsave(g4, filename = "Figures_v2/FigS1X_repStressGenes_4.pdf",
       height = 7.5, width = 12)

ga <- ggarrange(g1, g2, g3, g4, ncol = 1, nrow = 4)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.png",
       height = 30, width = 12)
ggsave(ga, filename = "Figures_v2/FigS1X_repStressGenes.pdf",
       height = 30, width = 12)



#####################################################################################
#################################### ssDRIP-Seq #####################################
#####################################################################################

cat("\n", timestamp2(), " Step 4: Generate ssDRIP-Seq figures\n", sep = "")

# Pre-process all samples
cat("\n", timestamp2(), " Preprocessing and alignment of ssDRIP-Seq reads...\n", sep = "")
runInfo <- read.csv("Data/DRIPSeq/Code/runInfo_2020.03.05_15.27.29.csv", stringsAsFactors = FALSE)
groupInfo <- read.table("Data/DRIPSeq/Code/groupFile_MultiOmics_ProgenitorStudy.txt", header = TRUE, stringsAsFactors = FALSE)
finalInfo <- merge(x = runInfo, y = groupInfo, by = "SampleName")
finalInfo$resultName <- paste0(finalInfo$SRAStudy, "_", 
                               finalInfo$Group, "_",
                               finalInfo$Condition, "_",
                               finalInfo$Replicate)
finalInfo <- finalInfo[finalInfo$Condition %in% c("S96", "Input"), c(1, 2, 34:40)]
runList <- finalInfo$Run
write.table(finalInfo, file = "Tables/ssDRIP_sampleInfo.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
dir.create("Data/DRIPSeq/tmp", showWarnings = FALSE)
system("export TMPDIR=Data/DRIPSeq/tmp")
dir.create("Data/DRIPSeq/Raw_Reads", showWarnings = FALSE)
dir.create("Data/DRIPSeq/Bam_Files/", showWarnings = FALSE)
for (i in 1:length(runList)) {
  currentRun <- runList[i]
  message("Processing SRA run ", currentRun)
  # Gather FASTQ files
  if (! file.exists(paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.fastq")) &
      ! file.exists(paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.qc.fastq")) &
      ! file.exists(paste0("Data/DRIPSeq/Bam_Files/", currentRun, "/", currentRun, ".bam"))) {
    cmd <- paste0("prefetch -O Data/DRIPSeq/tmp/ ", currentRun)
    downloaded <- FALSE
    while (! downloaded) {
      system(cmd)
      if (file.exists(paste0("Data/DRIPSeq/tmp/", currentRun, "/", currentRun, ".sra"))) {
        downloaded <- TRUE
      } else {
        downloaded <- FALSE
        warning("Failed to download SRA file. Trying again...")
      }
    }
    cmd <- paste0("parallel-fastq-dump --split-files -s ", paste0("Data/DRIPSeq/tmp/", currentRun, "/", currentRun, ".sra"), 
                  " -t ", maxCores, " -O Data/DRIPSeq/Raw_Reads/", currentRun, " --tmpdir Data/DRIPSeq/tmp &>/dev/null")
    system(cmd)
  } else {
    message("Fastq file found...")
  }
  # Run QC
  if (! file.exists(paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.qc.fastq")) &
      ! file.exists(paste0("Data/DRIPSeq/Bam_Files/", currentRun, "/", currentRun, ".bam"))) {
    cmd <- paste0("fastp -w ", maxCores,
                  " -i ", paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.fastq"), " -I ",
                  paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_2.fastq"), 
                  " -o ", paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.qc.fastq"),
                  " -O ", paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_2.qc.fastq"))
    system(cmd)
  } else {
    message("QC already performed...")
  }
  # Align reads
  if (! file.exists(paste0("Data/DRIPSeq/Bam_Files/", currentRun, "/", currentRun, ".bam"))) {
    dir.create(paste0("Data/DRIPSeq/Bam_Files/", currentRun), showWarnings = FALSE)
    cmd <- paste0("bwa mem -t ", maxCores, " Data/DRIPSeq/version0.6.0/genome.fa ",
                  paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_1.qc.fastq "),
                  paste0("Data/DRIPSeq/Raw_Reads/", currentRun, "/", currentRun, "_2.qc.fastq |",
                         " samtools view -b -@ ", maxCores, " - | samtools sort -@ ", maxCores,
                         " -o ", paste0("Data/DRIPSeq/Bam_Files/", currentRun, "/", currentRun, ".bam"),
                         " && samtools index ", paste0("Data/DRIPSeq/Bam_Files/", currentRun, "/", currentRun, ".bam")))
    system(cmd)
  } else {
    message("Alignment already performed...")
  }
}

# Split pos and negative alignments
cat("\n", timestamp2(), " Splitting alignments by first-in-pair strand orientation\n", sep = "")
availableBamsFiles <- list.files("Data/DRIPSeq/Bam_Files", recursive = TRUE, pattern = "\\.bam$", full.names = TRUE)
availableBamsFiles <- availableBamsFiles[file.size(availableBamsFiles) > 10000 &
                                           file.exists(paste0(availableBamsFiles, ".bai"))]
availableBams <- gsub(availableBamsFiles, pattern = "^Data/DRIPSeq/Bam_Files/.+/(.+)\\.bam$", replacement = "\\1")
names(availableBamsFiles) <- availableBams
dir.create("Data/DRIPSeq/bindingProfiles_singleStranded/", showWarnings = FALSE)
availableBamsFiles <- availableBamsFiles[names(availableBamsFiles) %in% paste0("SRR111853", c(67:90))]
for (i in 1:length(names(availableBamsFiles))) {
  availableBamNow <- availableBamsFiles[i]
  SRA <- names(availableBamNow)
  message(SRA)
  if (file.exists(paste0(gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg.bam")) &
      (file.size(paste0(gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg.bam")) > 10000)) {
    message("Already processed this bam file.")
    next
  }
  cmdF1 <- paste0("samtools view -f 99 -b ", availableBamNow, " -@ ", maxCores," -o ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_1.bam && ")
  cmdF2 <- paste0("samtools view -f 147 -b ", availableBamNow, " -@ ", maxCores," -o ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_2.bam && ")
  cmdF3 <- paste0("samtools merge -f -@ ", maxCores," ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos.bam ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_1.bam ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_2.bam && samtools index ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos.bam && ")
  cmdF <- paste0(cmdF1, cmdF2, cmdF3)
  cmdR1 <- paste0("samtools view -f 83 -b ", availableBamNow, " -@ ", maxCores," -o ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_1.bam && ")
  cmdR2 <- paste0("samtools view -f 163 -b ", availableBamNow, " -@ ", maxCores," -o ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_2.bam && ")
  cmdR3 <- paste0("samtools merge -f -@ ", maxCores," ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg.bam ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_1.bam ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_2.bam && samtools index ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg.bam && rm ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_1.bam && rm ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".pos_2.bam && rm ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_1.bam && rm ",
                  gsub(availableBamNow, pattern = "\\.bam$", replacement = ""), ".neg_2.bam")
  cmdR <- paste0(cmdR1, cmdR2, cmdR3)
  cmd <- paste0(cmdF, cmdR)
  system(cmd, wait = T)
}

# Call peaks
cat("\n", timestamp2(), " Peak calling...\n", sep = "")
for (i in 1:length(availableBamsFiles)) {
  availableBamNow <- availableBamsFiles[i]
  SRA <- names(availableBamNow)
  message(SRA)
  sampleNow <- finalInfo$SampleName[which(finalInfo$Run == SRA)]
  controlNow <- finalInfo$ControlSample[finalInfo$SampleName == sampleNow]
  if (is.na(controlNow)) {
    next
  }
  nameNow <- finalInfo$resultName[finalInfo$SampleName == sampleNow]
  dir.create(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow), showWarnings = FALSE)
  topBam <- availableBamNow
  ctrBam <- availableBamsFiles[finalInfo$Run[finalInfo$SampleName == controlNow]]
  topPos <- list.files(path = paste0("Data/DRIPSeq/Bam_Files/", names(topBam)), pattern = "\\.pos\\.bam$", full.names = TRUE)
  topNeg <- list.files(path = paste0("Data/DRIPSeq/Bam_Files/", names(topBam)), pattern = "\\.neg\\.bam$", full.names = TRUE)
  ctrPos <- list.files(path = paste0("Data/DRIPSeq/Bam_Files/", names(ctrBam)), pattern = "\\.pos\\.bam$", full.names = TRUE)
  ctrNeg <- list.files(path = paste0("Data/DRIPSeq/Bam_Files/", names(ctrBam)), pattern = "\\.neg\\.bam$", full.names = TRUE)
  cmd1 <- paste0("macs2 callpeak --nomodel --extsize 147 --broad -t ", topPos," -c ",
                 ctrPos, " -n ", paste0(nameNow, "_pos"),
                 " --outdir ", paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow), " && ")
  cmd2 <- paste0("macs2 callpeak --nomodel --extsize 147 --broad -t ", topNeg," -c ",
                 ctrNeg, " -n ", paste0(nameNow, "_neg"),
                 " --outdir ", paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow))
  cmd <- paste0(cmd1, cmd2)
  if (! file.exists(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow, "/", nameNow, "_neg_peaks.broadPeak"))) {
    system(cmd, wait = T)
  } else {
    message("Peaks already called for this sample.")
  }
}


# Process peaks into counts and run PCA
cat("\n", timestamp2(), " Assigning bam reads to peaks to construct peak count matrix...\n", sep = "")
finInfNow <- finalInfo
finInfNow$negTopBam <- sapply(finInfNow$Run, FUN = function(sampNow) {
  list.files(paste0("Data/DRIPSeq/Bam_Files/", sampNow), pattern = "\\.neg\\.bam$", full.names = TRUE)
})
finInfNow$posTopBam <- sapply(finInfNow$Run, FUN = function(sampNow) {
  list.files(paste0("Data/DRIPSeq/Bam_Files/", sampNow), pattern = "\\.pos\\.bam$", full.names = TRUE)
})
finInfNow2 <- finInfNow[! is.na(finInfNow$ControlType),]
finInfNow2$posCTRBam <- sapply(finInfNow2$Run, FUN = function(sampNow, finInfNow) {
  ctrSamp <- finInfNow2$ControlSample[finInfNow2$Run == sampNow]
  finInfNow$posTopBam[finInfNow$SampleName == ctrSamp]
}, finInfNow = finInfNow)
finInfNow2$negCTRBam <- sapply(finInfNow2$Run, FUN = function(sampNow, finInfNow) {
  ctrSamp <- finInfNow2$ControlSample[finInfNow2$Run == sampNow]
  finInfNow$negTopBam[finInfNow$SampleName == ctrSamp]
}, finInfNow = finInfNow)
finInfNow2$negPeaks <- sapply(finInfNow2$resultName, FUN = function(nameNow) {
  list.files(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow), pattern = "_neg_peaks.xls", full.names = TRUE)
})
finInfNow2$posPeaks <- sapply(finInfNow2$resultName, FUN = function(nameNow) {
  list.files(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", nameNow), pattern = "_pos_peaks.xls", full.names = TRUE)
})
# Build positive DBD object
posColData <- data.frame(
  "SampleID" = finInfNow2$resultName,
  "Tissue" = finInfNow2$Group,
  "Replicate" = finInfNow2$Replicate,
  "Batch" = ifelse(finInfNow2$Replicate %in% c(1, 2), 1, 2),
  "bamReads" = finInfNow2$posTopBam,
  "ControlID" = finInfNow2$ControlSample,
  "bamControl" = finInfNow2$posCTRBam,
  "Peaks" = finInfNow2$posPeaks,
  "PeakCaller" = "macs"
)
negColData <- data.frame(
  "SampleID" = finInfNow2$resultName,
  "Tissue" = finInfNow2$Group,
  "Replicate" = finInfNow2$Replicate,
  "Batch" = ifelse(finInfNow2$Replicate %in% c(1, 2), 1, 2),
  "bamReads" = finInfNow2$negTopBam,
  "ControlID" = finInfNow2$ControlSample,
  "bamControl" = finInfNow2$negCTRBam,
  "Peaks" = finInfNow2$posPeaks,
  "PeakCaller" = "macs"
)
if (! file.exists("Data/DRIPSeq/dbdPos_count.rda")) {
  dbdPos <- dba(sampleSheet = posColData)
  dbdPos <- dba.count(dbdPos)
  save(dbdPos, file = "Data/DRIPSeq/dbdPos_count.rda")
} else {
  cat("Already calculated counts! Loading and continuing...")
  load("Data/DRIPSeq/dbdPos_count.rda")
}

pdf("Figures_v2/DBA_PCA_Batch.pdf")
dba.plotPCA(dbdPos,attributes=c(DBA_TISSUE))
dev.off()

# Keep only batch one and try again with pos and neg
posColDataFinal <- posColData[posColData$Batch == 1 &
                                posColData$Tissue %in% c("MSCs", "iPSCs"),]
negColDataFinal <- negColData[negColData$Batch == 1 &
                                negColData$Tissue %in% c("MSCs", "iPSCs"),]

# Run DBA
if (! file.exists("Data/DRIPSeq/dbdPosFinal_count.rda")) {
  dbdPosFinal <- dba(sampleSheet = posColDataFinal)
  dbdPosFinal <- dba.count(dbdPosFinal)
  dbdPosFinal <- dba.contrast(dbdPosFinal, minMembers = 2)
  dbdPosFinal <- dba.analyze(dbdPosFinal)
  save(dbdPosFinal, file = "Data/DRIPSeq/dbdPosFinal_count.rda")
} else {
  cat("Already calculated counts! Loading and continuing...")
  load("Data/DRIPSeq/dbdPosFinal_count.rda")
}
if (! file.exists("Data/DRIPSeq/dbdnegFinal_count.rda")) {
  dbdnegFinal <- dba(sampleSheet = negColDataFinal)
  dbdnegFinal <- dba.count(dbdnegFinal)
  dbdnegFinal <- dba.contrast(dbdnegFinal, minMembers = 2)
  dbdnegFinal <- dba.analyze(dbdnegFinal)
  save(dbdnegFinal, file = "Data/DRIPSeq/dbdnegFinal_count.rda")
} else {
  cat("Already calculated counts! Loading and continuing...")
  load("Data/DRIPSeq/dbdnegFinal_count.rda")
}



# Keep only batch 1
countsNow <- dba.peakset(dbdPos, bRetrieve=TRUE)
sampsKeep <- which(posColData$Batch == 1)
countsNow <- countsNow[,sampsKeep]
posColData2 <- posColData[sampsKeep,]
countsNow <- as.matrix(mcols(countsNow))
pcsRes <- prcomp(t(countsNow))
all(colnames(countsNow) == posColData2$SampleID)
posColData2$PC1 <- pcsRes$x[,c(1)]
posColData2$PC2 <- pcsRes$x[,c(2)]
posColData2$Batch <- factor(posColData2$Batch)
perc1 <- round(pcsRes$sdev[1]/sum(pcsRes$sdev)*100)
perc2 <- round(pcsRes$sdev[2]/sum(pcsRes$sdev)*100)
gCols <- gg_color(6)
colMapDRIP <- c("iPSCs" = gCols[1],
                "hESCs" = gCols[2],
                "VECs" = gCols[3],
                "NSCs" = gCols[4],
                "MSCs" = gCols[5],
                "VSMCs" = gCols[6])
DimPlot2(data = posColData2, mapping = aes_string(x = "PC1", 
                                                 y = "PC2", 
                                                 color = "Tissue"),  
         xlab = paste0("PC1 (", perc1, "%)"), ylab = paste0("PC2 (", perc2, "%)"), 
         plotName = "DBA_PCA_Batch_Final", ptSize = 8, colorMap = colMapDRIP)

# Calculate number of peaks
cat("\n", timestamp2(), " Calculating numbers of peaks and peak overlaps...\n", sep = "")
if (! file.exists("Data/DRIPSeq/olapResList.rda")) {
  registerDoParallel(cores=6)
  samps <- list.dirs("Data/DRIPSeq/bindingProfiles_singleStranded/", 
                     recursive = F, full.names = FALSE)
  samps <- samps[grep(samps, pattern = "_1|_2")]
  peakList <- lapply(samps, FUN = function(sampNow) {
    # sampNow <- samps[1]
    # I had these backwards...
    posPeaks <- toGRanges(list.files(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", sampNow),
                                     pattern = "pos_peaks\\.xls", full.names = TRUE), format = "MACS2")
    colnames(mcols(posPeaks)) <- c(colnames(mcols(posPeaks))[-2], "name")
    strand(posPeaks) <- "-"
    negPeaks <- toGRanges(list.files(paste0("Data/DRIPSeq/bindingProfiles_singleStranded/", sampNow),
                                     pattern = "neg_peaks\\.xls", full.names = TRUE),  format = "MACS2")
    strand(negPeaks) <- "+"
    colnames(mcols(negPeaks)) <- c(colnames(mcols(negPeaks))[-2], "name")
    peaks <- c(posPeaks, negPeaks)
    peaksNow <- peaks[peaks$`-log10(qvalue)` > 3,]
    peaksNow
  })
  names(peakList) <- gsub(samps, pattern = "SRP.+_(.+)_S96_(.+)", replacement = "\\1_\\2")
  categories <- unique(gsub(samps, pattern = "SRP.+_(.+)_S96_(.+)", replacement = "\\1"))
  listList <- list()
  for (i in 1:length(categories)) {
    catNow <- categories[i]
    peakListCat <- peakList[grep(names(peakList), pattern = catNow)]
    peakListCat <- peakListCat[grep(names(peakListCat), pattern = "_1|_2")]
    listList[[i]] <- peakListCat
  }
  names(listList) <- categories
  runOLap <- function(listListNow) {
    olCat <- findOverlapsOfPeaks(listListNow, ignore.strand = FALSE)
    return(olCat)
  }
  olapResList <- foreach(k=1:6) %dopar% runOLap(listList[[k]])
  names(olapResList) <- categories
  save(olapResList, file = "Data/DRIPSeq/olapResList.rda")
} else {
  cat("\nOverlaps found! Loading and continuing...")
  load("Data/DRIPSeq/olapResList.rda")
}

mergePeakList <- list()
pltDF <- data.frame("SampleType" = rep(names(olapResList), each = 2),
                    "Replicate" = rep(c(1, 2), 6))
numPeaks <- c()
for (i in 1:length(olapResList)) {
  oLapNow <- olapResList[[i]]
  nameNow <- names(olapResList)[i]
  # makeVennDiagram(oLapNow)
  numPeaks <- c(numPeaks, sum(oLapNow$venn_cnt[,c(4)]), sum(oLapNow$venn_cnt[,c(5)]))
  namesNow <- as.character(oLapNow$peaksInMergedPeaks$name)
  peaksInMergedPeaks <- oLapNow$peaksInMergedPeaks[names(oLapNow$peaksInMergedPeaks) %in%
                                                     unlist(oLapNow$peaklist[[3]]$peakNames, 
                                                            use.names = FALSE),]
  peaksInMergedPeaksRed <- reduce(peaksInMergedPeaks, with.revmap = TRUE)
  peaksInMergedPeaksRed$fold_enrichment <- mean(extractList(peaksInMergedPeaks$fold_enrichment, peaksInMergedPeaksRed$revmap))
  peaksInMergedPeaksRed$qVal <- mean(extractList(peaksInMergedPeaks$`-log10(qvalue)`, peaksInMergedPeaksRed$revmap))
  mergePeakList[[nameNow]] <- peaksInMergedPeaksRed
}
pltDF$numPeaks <- numPeaks
pltDF <- pltDF[order(pltDF$numPeaks, decreasing = TRUE),]
pltDF$SampleType <- factor(pltDF$SampleType, levels = unique(pltDF$SampleType))
pal <-  colorRampPalette(brewer.pal(9,"Greys"))(length(unique(pltDF$SampleType)))
g1 <- ggplot(pltDF, mapping = aes_string(x = "SampleType",
                                         y = "numPeaks",
                                         fill = "SampleType")) +
  geom_boxplot() +  
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme_pubr(border = T, base_size = 22) +
  ylab("Number of DRIP peaks") +
  rremove("xlab") +
  rremove("legend") +
  scale_fill_manual(values = colMapDRIP) +
  stat_compare_means(comparisons = rev(list(c("iPSCs", "MSCs"),
                                            c("NSCs", "MSCs"),
                                            c("hESCs", "MSCs"))),
                     method = "t.test", label = "p.signif", size = 5, 
                     method.args = list("alternative" = "greater"))
ggsave(g1, filename = "Figures_v2/Num_DRIP_Peaks.pdf")

if (! file.exists("Data/DRIPSeq/oLapMerged2.rda")) {
  oLapMerged2 <- findOverlapsOfPeaks(mergePeakList[c(2, 3)], ignore.strand = FALSE)
  save(oLapMerged2, file = "Data/DRIPSeq/oLapMerged2.rda")
} else {
  load("Data/DRIPSeq/oLapMerged2.rda")
}

png(filename = "Figures_v2/DRIP_Peak_Overlap.png", units = "in",
    height = 5, width = 7, res = 300)
makeVennDiagram(oLapMerged2, cex = 1.5,
                fill = c(colMapDRIP[1], colMapDRIP[5]),
                ignore.strand = FALSE)
dev.off()
pdf(file = "Figures_v2/DRIP_Peak_Overlap.pdf", 
    height = 5, width = 7)
makeVennDiagram(oLapMerged2, cex = 1.5, 
                fill = c(colMapDRIP[1], colMapDRIP[5]),
                ignore.strand = FALSE)
dev.off()




suppressWarnings(rm("expression"))
suppressWarnings(rm("expressionNew"))
suppressWarnings(rm("vsd"))
suppressWarnings(rm("vsd"))
suppressWarnings(rm("countsNow"))
suppressWarnings(rm("phateNonEWS"))
suppressWarnings(rm("phateNow"))
suppressWarnings(rm("vsd2"))
suppressWarnings(rm("vsdTop"))
suppressWarnings(rm("vsdTop2"))
suppressWarnings(rm("vsdTopFinal"))
suppressWarnings(rm("vsdTopFinalNonEWS"))
suppressWarnings(rm("dbdPos"))
suppressWarnings(rm("olapResList"))
suppressWarnings(rm("oLapNow"))
suppressWarnings(rm("oLapMerged2"))
invisible(gc())


#####################################################################################
################################ Single Cell RNASeq #################################
#####################################################################################


# PHATE_1 signature in EWS datasets
if (! file.exists("Data/scRNASeq/BDInt.rda")) {
  # Add the Bishop group's single cell
  # -- CHLA9
  sm <- readMM("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368462_CHLA9_matrix.mtx.gz")
  samples <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368462_CHLA9_barcodes.tsv.gz", stringsAsFactors = FALSE)
  genesFinal <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368462_CHLA9_features.tsv.gz", stringsAsFactors = FALSE)
  genesFinalDups <- which(duplicated(genesFinal$V2))
  sm <- sm[-genesFinalDups,]
  genesFinal <- genesFinal$V2[-genesFinalDups]
  colnames(sm) <- samples$V1
  rownames(sm) <- genesFinal
  CHLA9 <- CreateSeuratObject(sm)
  CHLA9$Sample <- "CHLA9"
  # -- CHLA10
  sm <- readMM("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368463_CHLA10_matrix.mtx.gz")
  samples <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368463_CHLA10_barcodes.tsv.gz", stringsAsFactors = FALSE)
  genesFinal <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368463_CHLA10_features.tsv.gz", stringsAsFactors = FALSE)
  genesFinalDups <- which(duplicated(genesFinal$V2))
  sm <- sm[-genesFinalDups,]
  genesFinal <- genesFinal$V2[-genesFinalDups]
  colnames(sm) <- samples$V1
  rownames(sm) <- genesFinal
  CHLA10 <- CreateSeuratObject(sm)
  CHLA10$Sample <- "CHLA10"
  # -- TC71
  sm <- readMM("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368464_TC71_matrix.mtx.gz")
  samples <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368464_TC71_barcodes.tsv.gz", stringsAsFactors = FALSE)
  genesFinal <- read.table("Data/scRNASeq/CountMatrices/BISHOP_EWS_Cells/GSM4368464_TC71_features.tsv.gz", stringsAsFactors = FALSE)
  genesFinalDups <- which(duplicated(genesFinal$V2))
  sm <- sm[-genesFinalDups,]
  genesFinal <- genesFinal$V2[-genesFinalDups]
  colnames(sm) <- samples$V1
  rownames(sm) <- genesFinal
  TC71 <- CreateSeuratObject(sm)
  TC71$Sample <- "TC71"
  
  # -- Merge
  BISHOP <- merge(x = CHLA9, y = c(CHLA10, TC71))
  BISHOP$SRA <- "BISHOP"
  BISHOP$SRS <- BISHOP$Sample
  BISHOP$TissueSite <- "Ewing Sarcoma"
  BISHOP$chemistry <- "v3"
  BISHOP$Species <- "Homo sapiens"
  BISHOP$Counter <- "cellranger"
  # QC Figures
  BISHOP[["percent.mt"]] <- PercentageFeatureSet(BISHOP, pattern = "^MT-")
  pltData <- BISHOP@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Ai.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Aii.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "Sample", y = "nCount_RNA",
                                             color = "Sample")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Aiii.png", device = "png")
  # QC each
  BISHOP <- autoQCSRT(BISHOP, minNumCounts = 15000,
                      whichColumn = "SRS",
                      minNumGenes = 2500, maxPercentMT = 18 )
  # Show again after QC
  pltData <- BISHOP@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Aiv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Av.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures_v2/Fig2Avi.png", device = "png")
  # Normal clustering
  BISHOP <- NormalizeData(BISHOP)
  BISHOP <- FindVariableFeatures(BISHOP)
  BISHOP <- ScaleData(BISHOP)
  BISHOP <- RunPCA(BISHOP)
  BISHOP <- FindNeighbors(BISHOP, dims = 1:50)
  BISHOP <- FindClusters(BISHOP)
  BISHOP <- RunUMAP(BISHOP, dims = 1:50)
  BISHOP <- CellCycleScoring(BISHOP, 
                             s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
  pltData <- BISHOP@meta.data
  pltData$UMAP_1 <- BISHOP@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- BISHOP@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  DimPlot2(pltData, width = 12,
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Avii")
  # Needs batch correct 
  srtList <- SplitObject(BISHOP, split.by = "SRS")
  for (i in 1:length(srtList)) {
    srtList[[i]] <- NormalizeData(srtList[[i]], verbose = FALSE)
    srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = FALSE)
  }
  srtList <- FindIntegrationAnchors(object.list = srtList, dims = 1:30,
                                    anchor.features = 5000)
  BISHOPInt <- IntegrateData(anchorset = srtList, dims = 1:30)
  DefaultAssay(BISHOPInt) <- "integrated"
  BISHOPInt <- ScaleData(BISHOPInt)
  BISHOPInt <- RunPCA(BISHOPInt)
  BISHOPInt <- FindNeighbors(BISHOPInt)
  BISHOPInt <- FindClusters(BISHOPInt)
  BISHOPInt <- RunUMAP(BISHOPInt, dims = 1:50)
  # Separate plots
  pltData <- BISHOPInt@meta.data
  pltData$UMAP_1 <- BISHOPInt@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- BISHOPInt@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  pltData$Cluster <- pltData$seurat_clusters
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"), 
           plotName = "Fig2Aviii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Aix")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"), 
           plotName = "Fig2Ax")
  # # Check for expression of BACH1
  # FeaturePlot(BISHOP, features = "BACH1")
  # RidgePlot(BISHOP, features = "BACH1")
  
  
  # Add in the Delattre group's single cell
  filesNow <- list.files("Data/scRNASeq/CountMatrices/Delattre_EWS_PDX", pattern = "mtx", full.names = T)
  pdxIDs <- gsub(pattern = ".+_(.+)_matrix.+", x = filesNow, replacement = "\\1")
  names(filesNow) <- pdxIDs
  srtList <- lapply(names(filesNow), FUN = function(idNow, filesNow) {
    mtxNow <- filesNow[idNow]
    barcodesNow <- gsub(x = mtxNow, pattern = "matrix.mtx", replacement = "barcodes.tsv")
    genesNow <- gsub(x = mtxNow, pattern = "matrix.mtx", replacement = "genes.tsv")
    sm <- Matrix::readMM(mtxNow)
    genesNow <- read.table(genesNow)
    rownames(sm) <- gsub(genesNow$V2, pattern = "GRCh38_(.+)", replacement = "\\1")
    barcodesNow <- read.table(barcodesNow)
    colnames(sm) <- barcodesNow$V1
    srtTemp <- CreateSeuratObject(sm)
    srtTemp$SRA <- "Delattre"
    srtTemp$SRS <- idNow
    return(srtTemp)
  }, filesNow = filesNow)
  DELATTRE <- merge(x = srtList[[1]], y = c(srtList[c(2:length(srtList))]))
  DELATTRE$SRA <- "DELATTRE"
  DELATTRE$TissueSite <- "Ewing Sarcoma"
  DELATTRE$chemistry <- "v2"
  DELATTRE$Species <- "Homo sapiens"
  DELATTRE$Counter <- "cellranger"
  # QC Figures
  DELATTRE[["percent.mt"]] <- PercentageFeatureSet(DELATTRE, pattern = "^MT-")
  pltData <- DELATTRE@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Bi.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Bii.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Biii.png", device = "png")
  # QC each
  DELATTRE <- autoQCSRT(DELATTRE, minNumCounts = 4000,
                        maxNumCounts = 60000,
                        minNumGenes = 800, maxPercentMT = 20)
  # Show again after QC
  pltData <- DELATTRE@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Biv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Bv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 13.33,
         filename = "Figures_v2/Fig2Bvi.png", device = "png")
  
  # Normal clustering
  DELATTRE <- NormalizeData(DELATTRE)
  DELATTRE <- FindVariableFeatures(DELATTRE)
  DELATTRE <- ScaleData(DELATTRE)
  DELATTRE <- RunPCA(DELATTRE)
  DELATTRE <- FindNeighbors(DELATTRE)
  DELATTRE <- FindClusters(DELATTRE)
  DELATTRE <- RunUMAP(DELATTRE, dims = 1:50)
  DELATTRE <- CellCycleScoring(DELATTRE, 
                               s.features = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes)
  pltData <- DELATTRE@meta.data
  pltData$UMAP_1 <- DELATTRE@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- DELATTRE@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Bvii")
  # FeaturePlot(DELATTRE, features = "BACH1")
  # RidgePlot(DELATTRE, features = "BACH1")
  
  # Needs batch correct 
  srtList <- SplitObject(DELATTRE, split.by = "SRS")
  for (i in 1:length(srtList)) {
    srtList[[i]] <- NormalizeData(srtList[[i]], verbose = FALSE)
    srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = FALSE)
  }
  srtList <- FindIntegrationAnchors(object.list = srtList, dims = 1:30,
                                    anchor.features = 5000)
  DELATTREInt <- IntegrateData(anchorset = srtList, dims = 1:30)
  DefaultAssay(DELATTREInt) <- "integrated"
  DELATTREInt <- ScaleData(DELATTREInt)
  DELATTREInt <- RunPCA(DELATTREInt)
  DELATTREInt <- FindNeighbors(DELATTREInt)
  DELATTREInt <- FindClusters(DELATTREInt)
  DELATTREInt <- RunUMAP(DELATTREInt, dims = 1:30)
  
  # Separate plots
  pltData <- DELATTREInt@meta.data
  pltData$UMAP_1 <- DELATTREInt@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- DELATTREInt@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  pltData$Cluster <- pltData$seurat_clusters
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"), 
           plotName = "Fig2Bviii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Bix")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"), 
           plotName = "Fig2Bx")
  
  # Try combining BISHOP and Delattre
  BD <- merge(x = BISHOP, y = DELATTRE)
  BD <- NormalizeData(BD)
  BD <- ScaleData(BD, block.size = 10000)
  BD <- FindVariableFeatures(BD)
  BD <- RunPCA(BD)
  ElbowPlot(BD)
  BD <- FindNeighbors(BD, dims = 1:50)
  BD <- FindClusters(BD, resolution = .5)
  BD <- RunUMAP(BD, dims = 1:50)
  BD <- CellCycleScoring(BD, 
                         s.features = cc.genes$s.genes,
                         g2m.features = cc.genes$g2m.genes)
  pltData <- BD@meta.data
  pltData$UMAP_1 <- BD@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- BD@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Ci")
  # Needs batch correct 
  srtList <- SplitObject(BD, split.by = "SRS")
  for (i in 1:length(srtList)) {
    srtList[[i]] <- NormalizeData(srtList[[i]], verbose = FALSE)
    srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = FALSE)
  }
  srtList <- FindIntegrationAnchors(object.list = srtList, dims = 1:30,
                                    anchor.features = 5000)
  BDIntRaw <- IntegrateData(anchorset = srtList, dims = 1:30)
  DefaultAssay(BDIntRaw) <- "integrated"
  BDInt <- ScaleData(BDIntRaw)
  BDInt <- FindVariableFeatures(BDInt)
  BDInt <- RunPCA(BDInt)
  BDInt <- FindNeighbors(BDInt)
  BDInt <- FindClusters(BDInt)
  BDInt <- RunUMAP(BDInt, dims = 1:30)
  BDInt <- CellCycleScoring(BDInt, s.features = cc.genes$s.genes,
                            g2m.features = cc.genes$g2m.genes)
  pltData <- BDInt@meta.data
  pltData$UMAP_1 <- BDInt@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- BDInt@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  pltData$Cluster <- pltData$seurat_clusters
  pltData$Type <- "NONE"
  pltData$Type[pltData$SRA == "BISHOP"] <- "Cell line"
  pltData$Type[pltData$SRA == "DELATTRE"] <- "PDX"
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                color = "Cluster"), 
           plotName = "Fig2Cii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                color = "Sample"), 
           plotName = "Fig2Ciii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                color = "Phase"), 
           plotName = "Fig2Civ")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                color = "Type"), 
           plotName = "Fig2Cv")
  # Remove outlier cells
  pltData$Outlier <- F
  pltData$Outlier[pltData$UMAP_2 > 3.95] <- T
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                color = "Outlier"), 
           plotName = "Fig2Cvi")
  goodCells <- rownames(pltData)[! pltData$Outlier]
  BDInt <- subset(BDInt, cells = goodCells)
  # Reprocess
  BDInt <- ScaleData(BDInt)
  BDInt <- FindVariableFeatures(BDInt)
  BDInt <- RunPCA(BDInt)
  BDInt <- FindNeighbors(BDInt)
  BDInt <- FindClusters(BDInt)
  BDInt <- RunUMAP(BDInt, dims = 1:50)
  BDInt <- CellCycleScoring(BDInt, s.features = cc.genes$s.genes,
                            g2m.features = cc.genes$g2m.genes)
  save(BDInt, file = "Data/scRNASeq/BDInt.rda")
} else {
  cat("\nIntegrated Seurat object found! Loading and continuing...")
  load("Data/scRNASeq/BDInt.rda")
}

# Plot the integrated EWS dataset
DefaultAssay(BDInt) <- "integrated"
pltData <- BDInt@meta.data
pltData$UMAP_1 <- BDInt@reductions$umap@cell.embeddings[,c(1)]
pltData$UMAP_2 <- BDInt@reductions$umap@cell.embeddings[,c(2)]
pltData$Sample <- pltData$SRS
pltData$Cluster <- pltData$seurat_clusters
pltData$Type <- "NONE"
pltData$Type[pltData$SRA == "BISHOP"] <- "Cell line"
pltData$Type[pltData$SRA == "DELATTRE"] <- "PDX"
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                              color = "Type"), 
         plotName = "Fig2Cvii")
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                              color = "Cluster"), 
         plotName = "Fig2Cviii")
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                              color = "Sample"), 
         plotName = "Fig2Cix")
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                              color = "Phase"), 
         plotName = "Fig2Cx")

# Show cell cycle contribution in cell line and PDX
pltDataNoCell <- pltData[pltData$SRA == "DELATTRE",]
g1 <- ggplot(data = pltDataNoCell, mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                               color = "Phase")) +
  geom_point(size = .4) +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  scale_y_continuous(expand = c(0,0), limits = c(-5.25, 4.85))+
  scale_x_continuous(expand = c(0,0), limits = c(-7.8, 7.8))+
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  rremove("legend")
g1
ggsave(g1, filename = "Figures_v2/Phase_PDX.png", width = 8, height = 7.5)
pltDataCell <- pltData[pltData$SRA != "DELATTRE",]
g1 <- ggplot(data = pltDataCell, mapping = aes_string(x = "UMAP_1", y = "UMAP_2", 
                                                        color = "Phase")) +
  geom_point(size = .4) +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  scale_y_continuous(expand = c(0,0), limits = c(-5.25, 4.85))+
  scale_x_continuous(expand = c(0,0), limits = c(-7.8, 7.8))+
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  rremove("legend")
g1
ggsave(g1, filename = "Figures_v2/Phase_Cells.png", width = 8, height = 7.5)




# Cluster heatmap
tableDF <- as.data.frame(table(pltData$SRS, pltData$Cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
ph <- pheatmap(t(resMat), cluster_rows = T, silent = TRUE,
         cluster_cols = F, labels_row = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 18)
plt <- ggplotify::as.ggplot(ph[[4]])
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures_v2/Fig2Cxi.png", device = "png")
tableDF <- as.data.frame(table(pltData$Phase, pltData$Cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
ph <- pheatmap(t(resMat), cluster_rows = T, silent = TRUE,
               cluster_cols = F, labels_row = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 18)
plt <- ggplotify::as.ggplot(ph[[4]])
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures_v2/Fig2Cxii.png", device = "png")

if (! file.exists("Data/scRNASeq/BDIntMarkers.rda")) {
  # Find cluster markers
  plan("multiprocess", workers = length(unique(BDInt@meta.data$seurat_clusters)))
  DefaultAssay(BDInt) <- "integrated"
  BDIntMarkers <- FindAllMarkers(BDInt, logfc.threshold = .125)
  save(BDIntMarkers,  file = "Data/scRNASeq/BDIntMarkers.rda")
  plan("sequential")
  doClusterMarkerPlot(srt = BDInt, topN = 10, intClust = 8,
                      TERM2GENE = TERM2GENE_CGP,
                      srtMarkers = BDIntMarkers, 
                      outName = "Fig2Cxiii")
  doClusterMarkerPlot(srt = BDInt, TERM2GENE = TERM2GENE_EWS,
                      srtMarkers = BDIntMarkers, intClust = 10,
                      outName = "Fig2Cxiv")
} else {
  load("Data/scRNASeq/BDIntMarkers.rda")  
}

write.table(BDIntMarkers, file = "Tables/clusterMarkers_BDInt.tsv", 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Plot PHATE_1 effect from bulk RNA-Seq
DefaultAssay(BDInt) <- "RNA"
BDInt <- NormalizeData(BDInt)  
BDInt <- ScaleData(BDInt, do.center = F, 
                   block.size = 10000)
load("Data/bulkRNASeq/trendyRes_NormalPHATE_WithoutEWSTrendy.rda")
rankDF <- resListWithoutEWS$PHATE_1 
countsNow <- BDInt@assays$RNA@scale.data
ranks <- rankDF$value
names(ranks) <- rankDF$geneName
ranks <- ranks[names(ranks) %in% rownames(countsNow)]
countsNow <- countsNow[names(ranks),]
countsNowPHATE <- countsNow * ranks
phate_1_score <- apply(countsNowPHATE, MARGIN = 2, FUN = function(colNow) {
  median(colNow[abs(colNow) > 0])
})

BDInt$PHATE_1_SCORE <- phate_1_score
pltData <- BDInt@meta.data
keepSamps <- c()
for (i in 1:length(unique(pltData$seurat_clusters))) {
  datNow <- pltData[pltData$seurat_clusters == unique(pltData$seurat_clusters)[i],]
  p1Outs <- boxplot.stats(datNow$PHATE_1_UP_SCORE)$out
  keepSamps <- c(keepSamps, rownames(datNow)[! datNow$PHATE_1_UP_SCORE %in% p1Outs])
}
pltDataNow <- pltData[rownames(pltData) %in% keepSamps,]
pltData$Cluster <- paste0("Cluster_", pltData$seurat_clusters)
pltData <- pltData[order(pltData$seurat_clusters),]
pltData$Cluster <- factor(pltData$Cluster, levels = unique(pltData$Cluster))
g1 <- ggplot(pltData, mapping = aes_string(x = "PHATE_1_SCORE",
                                           y = "Cluster",
                                           fill = "Cluster")) +
  geom_density_ridges(scale = 2) +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  rremove("legend") + rremove("ylab") + 
  xlab("Developmental Context PHATE_1 score")
ggsave(g1, filename = "Figures_v2/Fig2Cxiv.png", height = 7.5, width = 9)

# EWS Markers highlight EWS cells -- some other tissues have higher expression of these
pltData <- BDInt@meta.data
pltData$UMAP_1 <- BDInt@reductions$umap@cell.embeddings[,c(1)]
pltData$UMAP_2 <- BDInt@reductions$umap@cell.embeddings[,c(2)]
pltData$Cluster <- paste0("Cluster_", pltData$seurat_clusters)
meds <- sapply(unique(pltData$Cluster), FUN = function(Cluster) {
  median(pltData$PHATE_1_SCORE[pltData$Cluster == Cluster])
})
names(meds) <- unique(pltData$Cluster)
meds <- meds[order(meds)]
pltData <- pltData[order(match(pltData$Cluster, names(meds))),]
pltData$Cluster <- factor(pltData$Cluster, levels = unique(pltData$Cluster))
g1 <- ggplot(pltData, 
             mapping = aes_string(x = "Cluster", fill = "Cluster",
                                  y = "PHATE_1_SCORE")) +
  geom_boxplot() +
  theme_pubr(border = T, base_size = 22) +
  rremove("legend") + 
  rremove("ylab") +
  ylab("PHATE_1 Score") +
  ggpubr::rotate() +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis(discrete = T)
ggsave(g1, filename = "Figures_v2/PHATE_1_SCORE_scRNASEQ_Clusters.png",
       height = 7.5, width = 9)

pltDataSave <- pltData[,c(5, 6, 2, 3, 11, 16, 22, 19, 20, 21)]
colnames(pltDataSave)[c(1:2)] <- c("Study", "Sample")
write.table(pltDataSave, file = "Tables/scRNASeq_BDInt_FinalInfo.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

doClusterMarkerPlot(srt = BDInt, intClust = c(8),
                     TERM2GENE = TERM2GENE,
                     srtMarkers = BDIntMarkers, 
                     outName = "IntClust")
BDIntMarkers
pltDat <- BDInt@meta.data
pltDat$UMAP_1 <- BDInt@reductions$umap@cell.embeddings[,c(1)]
pltDat$UMAP_2 <- BDInt@reductions$umap@cell.embeddings[,c(2)]
clusterNow <- 8
posSelect <- BDIntMarkers$gene[BDIntMarkers$cluster == 8 &
                                 BDIntMarkers$p_val_adj < .05 &
                                 BDIntMarkers$avg_logFC > 0]
negSelect <- BDIntMarkers$gene[BDIntMarkers$cluster == 8 &
                                 BDIntMarkers$p_val_adj < .05 &
                                 BDIntMarkers$avg_logFC < 0]
PathRes <- doPathEnrichPlot(genesUp = posSelect, genesDn = negSelect, TERM2GENE = TERM2GENE, returnData = TRUE)
eresUp <- PathRes$eresUp
eresUp$Group <- "Over-expressed"
eresDn <- PathRes$eresDn
eresDn$Group <- "Under-expressed"
eres <- rbind(eresUp[,c(-1)], eresDn[c(-1)])
eres <- eres[,c(9, 1:8)]
write.table(eres, file = "Tables/Cluster_8_Marker_Enrichment.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
pathResPlot <- PathRes$plot
ggsave(pathResPlot, height = 6, width = 9, filename = "Figures_v2/Cluster_8_pathEnrichplot.png")


# ## Do BACH1 markers ##
# bachMarkers <- c("MMP1", "CXCR4", "MMP9", "MMP13", "HMGA2", "VIM", "BCL2",
#                  "BCL2L1", "HMOX1", "CCND1", "GCLM", "FTH1", "PSAP", "SOCS2",
#                  "CALM1", "COPS6", "CDK6", "MAFG", "EWSR1", "LRRC8D", "BCL2L11",
#                  "SQSTM1", "RHBDD3", "RNFRSF1A", "AFG3L1", "ZNF3", "PEBP1",
#                  "ME1", "ALDOA", "TKT", "SLC7A11", "SLC25A10", "MAPT", "SLC48A1",
#                  "CLSTN1", "IFNB1", "IRF9", "IFITM1", "MX1", "OAS1", "OAS2",
#                  "IL6ST", "STEAP3")
# # Get counts
# countsNorm <- BDInt@assays$RNA@scale.data
# pltData <- BDInt@meta.data
# pltData$UMAP_1 <- BDInt@reductions$umap@cell.embeddings[,c(1)]
# pltData$UMAP_2 <- BDInt@reductions$umap@cell.embeddings[,c(2)]
# pltData$Cluster <- paste0("Cluster_", pltData$seurat_clusters)
# 
# g1 <- featurePlot2(pltData = pltData, countsNorm = countsNorm,
#                    featureName = "BACH1-related", plotName = "BACH1_Markers",
#                    features = bachMarkers, x = "UMAP_1", y = "UMAP_2")
# g1
# ga <- FeaturePlot(BDInt, features = bachMarkers, ncol = 4)
# ggsave(ga, filename = "BACH1_related_featurePlot.png", limitsize = F,
#        height = 70, width = 25)
# 
# 
# countsSmall <- countsNorm[rownames(countsNorm) %in% bachMarkers,]
# dim(countsSmall)
# countsSmallPR <- prcomp(countsSmall)
# newPlt <- data.frame("PC1" = countsSmallPR$x[,c(1)],
#                      "PC2" = countsSmallPR$x[,c(2)], stringsAsFactors = FALSE,
#                      "geneName" = rownames(countsSmallPR$x))
# community <- FindNeighbors(countsSmallPR$x)
# clusts <- FindClusters(community$snn, resolution = 1.04)
# newPlt$Cluster <- clusts$res.1.04
# ggscatter(newPlt, x = "PC1", label = "geneName", legend = "none", repel = TRUE,
#           y = "PC2", color = "Cluster")
# 
# hypoxicBACH1Genes <- newPlt$geneName[newPlt$Cluster == 0]
# cyclingBACH1Genes <- newPlt$geneName[newPlt$Cluster == 1]
# ubiquioutBACH1Genes <- newPlt$geneName[newPlt$Cluster == 2]
# g1 <- featurePlot2(pltData = pltData, countsNorm = countsNorm,
#                    featureName = "BACH1-related (Hypoxia)", plotName = "BACH1_Markers_Hypoxia",
#                    features = hypoxicBACH1Genes, x = "UMAP_1", y = "UMAP_2")
# g1
# g1 <- featurePlot2(pltData = pltData, countsNorm = countsNorm,
#                    featureName = "BACH1-related (Cycling)", plotName = "BACH1_Markers_Cycling",
#                    features = cyclingBACH1Genes, x = "UMAP_1", y = "UMAP_2")
# g1
# 
# 
# g1 <- featurePlot2(pltData = pltData, countsNorm = countsNorm,
#                    featureName = "BACH1-related (Ubiquitous)", plotName = "BACH1_Markers_Ubiquitous",
#                    features = ubiquioutBACH1Genes, x = "UMAP_1", y = "UMAP_2")
# g1


cat(timestamp2(), " DONE!")

