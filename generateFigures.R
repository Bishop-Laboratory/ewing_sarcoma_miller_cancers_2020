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

# TODO: Script for generating the figures in DOI:??

# For reproducibility:
## TODO: 1. Clone github repo at:
## 2. Install conda https://docs.conda.io/en/latest/miniconda.html
## TODO: 3. cd NAME_OF_GIT_REPO/
## 4. conda env create -f environment.yml
## 5. conda activate ewsPaperEnv
## 6. (Rscript Code/generateFigures.R) |& tee Code/generateFigures_logFile.txt

#####################################################################################
#####################################################################################

cat("\nRunning figure generation script for DOI:???\n")
source("helpers.R") 

#####################################################################################
############################ Preliminary: Load libraries ############################ 
#####################################################################################

cat("\n", timestamp2(), " Step 1: Loading libraries\n", sep = "")
suppressMessages(library(ggpubr))
suppressMessages(library(rhdf5))
suppressMessages(library(DESeq2))
suppressMessages(library(uwot))
suppressMessages(library(biomaRt))
suppressMessages(library(phateR))
suppressMessages(library(tidyr))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(grid))
suppressMessages(library(ggplotify))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(wordcloud))
suppressMessages(library(plotly))
suppressMessages(library(Trendy))
suppressMessages(library(clusterProfiler))
suppressMessages(library(EBSeq))
suppressMessages(library(ggrepel))
suppressMessages(library(VennDiagram))
suppressMessages(library(Seurat))
suppressMessages(library(ontologyIndex))
suppressMessages(library(jsonlite))
# Reproducible...
set.seed(42)

#####################################################################################
########################### Preliminary: Download Datasets ########################## 
#####################################################################################

cat("\n", timestamp2(), " Step 2: Downloading any missing datasets\n", sep = "")
if (! file.exists("Data/bulkRNASeq/human_matrix.h5")) {
  cat("\nARCHS4 expression data: \n")
  download.file(url = "https://mssm-seq-matrix.s3.amazonaws.com/human_matrix.h5", 
                destfile = "Data/bulkRNASeq/human_matrix.h5")
}
if (! file.exists("Data/bulkRNASeq/bto.obo")) {
  cat("\nBTO ontology file: \n")
  download.file(url = "https://raw.githubusercontent.com/BRENDA-Enzymes/BTO/master/bto.obo", 
                destfile = "Data/bulkRNASeq/bto.obo")
}

#####################################################################################
#################################### Bulk RNASeq ####################################
#####################################################################################

cat("\n", timestamp2(), " Step 3: Generate Bulk RNA-Seq figures\n", sep = "")
cat("\nFiltering unwanted samples and categorizing expression metadata...\n")

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
  
  # Filter out single cell datasets
  patternSingle <- c("single cell", "single-cell", "smart seq", "in-drop",
                     "cel-seq", "10X", "scRNA seq", "smartseq", "CELseq",
                     "smart-seq", "indrop", "drop-seq", "drop seq", "single nucleus",
                     "single-nucleus", "snRNA-Seq", "snRNASeq",
                     "fluidigm", "scRNASeq", "scRNA-Seq", "chromium")
  colData <- colData[unique(grep(pattern = paste(patternSingle, collapse="|"), 
                                 invert = T, perl = T,
                                 x = colData$extractProtocol, ignore.case = T)),]
  
  # Categorize cancer samples
  cancer = c("hela", "hek", "k562", "reh", "jurkat", "leukemi", "293", "bewo",
             "kras", "mcf", "lncap", "bjab", "gbm", "rko", "ramos", "mel888",
             "vcap", "saos2", "vapc", "nalm6", "set2", "tov21",
             "cancer", "carcin", "sarcom", "metasta", "tumor",
             "[a-zA-Z]oma ", "[a-zA-Z]oma$")
  cl <- ontologyIndex::get_OBO("Data/bulkRNASeq/bto.obo")
  cancerLines <- cl$name
  cancer <- c(cancer,
              as.character(cancerLines[grep(x = cancerLines, ignore.case = T, perl = T,
                                            pattern = paste0(cancer, collapse = "|"))]))
  colData$disease <- "Normal"
  colData$disease[grep(x = colData$tissue, pattern = paste0(cancer, collapse = "|"),
                       ignore.case = T, perl = T)] <- "Tumor"
  load("Data/bulkRNASeq/poorlyDescribedTumorSamples.rda")
  colData$disease[colData$samples %in% missingTumorSamps] <- "Tumor"
  
  # Categorize remaining samples
  tissueDict <- read_json("Data/bulkRNASeq/tissueDictionary.json")
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
  
  # Filter for read alignments > 10 million reads
  colDataFinal <- colDataFinal[colDataFinal$readsAligned > 10E6,]
  save(colDataFinal, file = "Data/bulkRNASeq/colDataFinal.rda")
} else {
  cat("\ncolData File Found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/colDataFinal.rda")
}

cat("\nLoading & filtering expression data...\n")
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

cat("\nApplying Variance Stabilizing Transform to dataset...\n")
if (! file.exists("Data/bulkRNASeq/fullVSTCounts.rda")) {
  # Make data Homoscedastic with VST 
  vsd <- vst(expression, nsub = 1000)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/bulkRNASeq/fullVSTCounts.rda")
} else if (file.exists("Data/bulkRNASeq/fullUMAPData.rda")) {
  cat("\nUMAP results found... Skipping this step...\n")
} else {
  cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullVSTCounts.rda")
}

cat("\nCalculating UMAP...\n")
if (! file.exists("Data/bulkRNASeq/fullUMAPData.rda")) {
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
    pcData <- pcData$x
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
  cat("\nUMAP data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullUMAPData.rda")
}

# Now without the EWS genes
cat("\nRecalculating VST with EWS geneset subtracted...\n")
if (! file.exists("Data/bulkRNASeq/nonEWSVSTCounts.rda")) {
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  load("Data/ewsGenes.rda")
  # Filter out the left side of the umap
  keepSamps <- pltData$samples
  nonEWSGenes <- rownames(expression)[! rownames(expression) %in% ewsGenes]
  expression <- expression[nonEWSGenes,keepSamps]
  # Make data Homoscedastic with VST
  vsd <- vst(expression, nsub = 1000)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/bulkRNASeq/nonEWSVSTCounts.rda")
} else if (file.exists("Data/bulkRNASeq/nonEWSUMAPData.rda")) {
  cat("\nUMAP results found... Skipping this step...\n")
} else {
  cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/nonEWSVSTCounts.rda")
}

cat("\nRecalculating UMAP with EWS geneset subtracted...\n")
if (! file.exists("Data/bulkRNASeq/nonEWSUMAPData.rda")) {
  load("Data/bulkRNASeq/colDataFinal.rda")
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
  if (! file.exists( "Data/bulkRNASeq/nonEWSPCData.rda")) {
    pcDataFull <- prcomp(t(vsdHVG))
    pcData <- pcData$x
    pcData <- pcData[,c(1:100)]
    save(pcData, file = "Data/bulkRNASeq/nonEWSPCData.rda")
  } else {
    load("Data/bulkRNASeq/nonEWSPCData.rda")
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
  save(pltData, file = "Data/bulkRNASeq/nonEWSUMAPData.rda")
} else {
  cat("\nUMAP data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/nonEWSUMAPData.rda")
}

# cat("\nCalculating cell type marker genes...\n")
# if (! file.exists("Data/bulkRNASeq/cellMarkerGenes.rda")) {
#   # Load full UMAP/VST data
#   load("Data/bulkRNASeq/fullUMAPData.rda")
#   load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
#   all(pltData$samples == colnames(expression))
#   srt <- CreateSeuratObject(expression, meta.data = pltData)
#   srt <- NormalizeData(srt)
#   srt <- ScaleData(srt)
#   srt$TissueType <- pltData$TissueType
#   Idents(srt) <- srt$TissueType
#   srtMarkers <- FindAllMarkers(srt, only.pos = TRUE, random.seed = 42)
#   save(srtMarkers, file = "Data/bulkRNASeq/cellMarkerGenes.rda")
#   save(srt, file = "Data/bulkRNASeq/srt.rda")
# } else {
#   cat("\nCell type marker genes found! Loading and continuing...\n")
#   load("Data/bulkRNASeq/cellMarkerGenes.rda")
# }


cat("\nGenerating bulk RNA-Seq figures...\n")
# Get cell type markers from panglao DB
load("Data/ewsCellMarkers.rda")
cellTypeMarkers <- read.table('Data/PanglaoDB_markers_07_Feb_2020.tsv', 
                              sep = "\t", header = T, stringsAsFactors = F)
cellTypeMarkers <- cellTypeMarkers[grep(cellTypeMarkers$species, 
                                        pattern = "Hs", perl = T),]
immunoMarkers <- cellTypeMarkers$official.gene.symbol[grep(cellTypeMarkers$organ, pattern = "Immune")]
endoMarkers <- cellTypeMarkers$official.gene.symbol[grep(cellTypeMarkers$cell.type, pattern = "Endothe")]
adiposeMarkers <- cellTypeMarkers$official.gene.symbol[grep(cellTypeMarkers$cell.type, pattern = "Adipo")]
ewingMarkers <- c(ewsCellMarkers)

# Load full UMAP/VST data
load("Data/bulkRNASeq/fullUMAPData.rda")
load("Data/bulkRNASeq/fullVSTCounts.rda")
immunoMarkers <- immunoMarkers[immunoMarkers %in% rownames(vsd)]
pltData$immuneScore <- scale(colMedians(vsd[immunoMarkers,]), center = F, scale = T)
endoMarkers <- endoMarkers[endoMarkers %in% rownames(vsd)]
adiposeMarkers <- adiposeMarkers[adiposeMarkers %in% rownames(vsd)]
pltData$adiposeScore <- scale(colMedians(vsd[adiposeMarkers,]), center = F, scale = T)
pltData$endoScore <- scale(colMedians(vsd[endoMarkers,]), center = F, scale = T)
ewingMarkers <- ewingMarkers[ewingMarkers %in% rownames(vsd)]
pltData$ewsScore <- scale(colMedians(vsd[ewingMarkers,]), center = F, scale = T)
pltData$TissueType[pltData$TissueType %in% c("HSCs",
                                             "MSCs",
                                             "neural crest cells")] <- "stem-like"
pltData$Cluster <- pltData$cluster
# # Identify probable latent tumor or oncogenic samples not caught by classifier
# # Only 40 samples identified, not necessary to re-run the full pipeline... 
# # Will simply exclude from downstream...
# pltData$latentTumor <- F
# pltData$latentTumor[pltData$UMAP_1 < 4 & pltData$UMAP_1 > -1 &
#                       pltData$UMAP_2 < 5 & pltData$UMAP_2 > -1.2] <- T
# latentTumorCells <- pltData$samples[(pltData$ewsScore > 1.4 | 
#                                        pltData$latentTumor) &
#                                       pltData$Tissue != "Ewing sarcoma"]
# pltData <- pltData[! pltData$samples %in% latentTumorCells,]
# vsd <- vsd[, colnames(vsd) %in% pltData$samples]
# Start plotting
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Ai.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Ai_legend.png", device = "png", height = 7.5, width = 12)
pltData$Tissue <- stringr::str_to_sentence(pltData$TissueType)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Aii.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Aii_legend.png", device = "png", height = 7.5, width = 12)
pltData$Tissue[! pltData$Tissue %in% c("Ewing sarcoma", "Stem-like")] <- "Mixed/Other"
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Aiii.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Aiii_legend.png", device = "png", height = 7.5, width = 12)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "immuneScore")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + labs(title = "Immune score") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "Greens"))(255))
ggsave(g1, filename = "Figures/Fig1Ei.png", device = "png", height = 8, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "endoScore")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + labs(title = "Endothelial score") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "Blues"))(255))
ggsave(g1, filename = "Figures/Fig1Eii.png", device = "png", height = 8, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "ewsScore")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + labs(title = "Ewing score") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "Reds"))(255))
ggsave(g1, filename = "Figures/Fig1Eiii.png", device = "png", height = 8, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "adiposeScore")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  rremove("legend") + labs(title = "Adipose score") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_gradientn(colors = colorRampPalette(brewer.pal(9, "Oranges"))(255))
ggsave(g1, filename = "Figures/Fig1Eiv.png", device = "png", height = 8, width = 8)


# zoom in on Clusters to merge
pltDataSmall <- pltData[pltData$cluster %in% c(2, 4, 14),]
pltDataSmall$Cluster <- "Undefined"
pltDataSmall$Cluster[pltDataSmall$cluster == 14] <- "14 (Ewing sarcoma)"
pltDataSmall$Cluster[pltDataSmall$cluster == 4] <- "4 (MSCs & Fibroblasts)"
pltDataSmall$Cluster[pltDataSmall$cluster == 2] <- "2 (Pluripotent & Neuronal progenitors)"

g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_color_manual(values = c("#619CFF", "#DA8F00","#A9A400")) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Bi.png", device = "png")
g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_color_manual(values = c("#619CFF", "#DA8F00","#A9A400")) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Bi_legend.png", device = "png")

pltData$showClust <- "No"
pltData$showClust[pltData$cluster %in% c(2, 4, 14)] <- "Yes"
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "showClust")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(linetype = 3, colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = c("#c4bfbe", "#e83c1a")) +
  rremove("legend") + rremove("axis") +
  rremove("ticks") + rremove("xy.text") +
  rremove("xylab")
ggsave(g1, height = 7.32, width = 7.32,
       filename = "Figures/Fig1Bi_inset.png", device = "png")

load("Data/bulkRNASeq/nonEWSUMAPData.rda")
pltData$TissueType[pltData$TissueType %in% c("HSCs",
                                             "MSCs",
                                             "neural crest cells")] <- "stem-like"
pltData$Cluster <- pltData$cluster
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Aiv.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Aiv_legend.png", device = "png", height = 7.5, width = 12)

pltData$Tissue <- stringr::str_to_sentence(pltData$TissueType)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Av.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Av_legend.png", device = "png", height = 7.5, width = 12)

pltData$Tissue[! pltData$Tissue %in% c("Ewing sarcoma", "Stem-like")] <- "Mixed/Other"
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, filename = "Figures/Fig1Avi.png", device = "png", height = 7.5, width = 8)
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, filename = "Figures/Fig1Avi_legend.png", device = "png", height = 7.5, width = 12)


# Classify
tissueDictNow <- jsonlite::read_json("Data/bulkRNASeq/tissueDictionary.json")
pltDataToClass <- unique(pltData[,c(1:3)])
pltDataToClass <- categorizeMetaData(pltDataToClass,
                                     cols = "tissue", 
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
catClass$Tissue[is.na(catClass$Tissue)] <- "Undefined"
pltData$TissueFinal <- catClass$Tissue
pltDataSmall <- pltData[pltData$cluster %in% c(1),]
pltDataSmall$Tissue <- "Undefined"
pltDataSmall$Tissue[! pltDataSmall$TissueFinal %in% c("stem-like",
                                                      "brain", "hESCs", 
                                                      "Neural progenitor/stem cells",
                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                      "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltDataSmall$Tissue[ pltDataSmall$TissueFinal %in% c("stem-like",
                                                     "brain", "hESCs", 
                                                     "Neural progenitor/stem cells",
                                                     "fibroblasts", "iPSCs",  "MSCs",
                                                     "ewing sarcoma", "neural crest cells")] <- pltDataSmall$TissueFinal[ pltDataSmall$TissueFinal %in% c("stem-like",
                                                                                                                                                          "brain", "hESCs", 
                                                                                                                                                          "Neural progenitor/stem cells",
                                                                                                                                                          "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                                          "ewing sarcoma", "neural crest cells")]
pltDataSmall$Tissue[pltDataSmall$Tissue %in% c("neural crest cells", "brain", "Neural progenitor/stem cells")] <- "Neural progenitors & brain"
pltDataSmall$Tissue[! pltDataSmall$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltDataSmall$Tissue[! pltDataSmall$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
pltDataSmall$Cluster <- "1 (Formerly 2, 4, 14)"
g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Cluster")) + 
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_color_manual(values = c("#EB8335", "#619CFF")) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Bii.png", device = "png")
g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Cluster")) + 
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  scale_color_manual(values = c("#EB8335", "#619CFF")) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Bii_legend.png", device = "png")
# Display tissue types found within new merged cluster
g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Biii.png", device = "png")

g1 <- ggplot(pltDataSmall, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                                color = "Tissue")) +
  geom_point(size = .5) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Biii_legend.png", device = "png")

# Over-representation heatmaps
load("Data/bulkRNASeq/fullUMAPData.rda")
classMap <- toClass[,c(1, 2, 3, 40)]
pltData <- pltData[,c(-13)]
catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
catClass$Tissue[is.na(catClass$Tissue)] <- "Undefined"
catClass <- catClass[catClass$Tissue != "Non-Ewing tumor",]
tableDF <- as.data.frame(table(catClass$Tissue, catClass$cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resDF <- resDF[rowSums(resDF) > 100,]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = T, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")] <- stringr::str_to_sentence(colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")])

pheatmap(resMat, cluster_rows = F, cluster_cols = T, labels_col = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)

ap <- grid.grab()

plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig1Ci_legend.png", device = "png")
load("Data/bulkRNASeq/nonEWSUMAPData.rda")
classMap <- toClass[,c(1, 2, 3, 40)]
pltData <- pltData[,c(-13)]
catClass <- merge(x = pltData, y = classMap, by = c("samples", "series", "tissue"), all.x = T)
catClass <- catClass[order(match(catClass$samples, pltData$samples)),]
catClass$Tissue[is.na(catClass$Tissue)] <- "Undefined"
catClass <- catClass[catClass$Tissue != "Non-Ewing tumor",]
tableDF <- as.data.frame(table(catClass$Tissue, catClass$cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resDF <- resDF[rowSums(resDF) > 100,]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = T, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")] <- stringr::str_to_sentence(colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")])

pheatmap(resMat, cluster_rows = F, cluster_cols = T, labels_col = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)
ap <- grid.grab()
plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig1Cii_legend.png", device = "png")

# Just cluster 1 now!
clusterOneColData <- catClass[catClass$cluster == 1,]
tableDF <- as.data.frame(table(clusterOneColData$Tissue))
tablDF2 <- as.data.frame(table(catClass$Tissue))
mergeTable <- merge(x = tableDF, y = tablDF2, by = "Var1")
# Only keep tissue cats in which at least 20% of total pop in cluster
goodCats <- as.character(mergeTable$Var1[(mergeTable$Freq.x/mergeTable$Freq.y) >= .2])
clusterOneColData <- catClass[catClass$cluster == 1 & catClass$Tissue %in% goodCats &
                                ! catClass$Tissue %in% "Non-Ewing tumor",]


# Cluster 1 with the EWS genes
cat("\nRe-calculating VST for samples from cluster 1...\n")
if (! file.exists("Data/bulkRNASeq/fullVSTCounts_clusterOne.rda")) {
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  keepSamps <- clusterOneColData$samples
  expression <- expression[,keepSamps]
  # Make data Homoscedastic with VST
  vsd <- vst(expression, nsub = 1000)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/bulkRNASeq/fullVSTCounts_clusterOne.rda")
} else if (file.exists("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")) {
  cat("\nUMAP results found... Skipping this step...\n")
} else {
  cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullVSTCounts_clusterOne.rda")
}

cat("\nCalculating UMAP with EWS genes...\n")
if (! file.exists("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")) {
  # Calculate highly-variable genes (HVGs)
  rv <- matrixStats::rowVars(vsd)
  names(rv) <- rownames(vsd)
  rv <- rv[order(rv, decreasing = T)]
  hvgs <- rv[c(1:10000)]
  vsdHVG <- vsd[names(hvgs),]
  # dim(vsdHVG)
  rm(vsd)
  gc()
  colDataNow <- clusterOneColData[order(match(clusterOneColData$samples, colnames(vsdHVG))),]
  if (! all(colDataNow$samples == colnames(vsdHVG))) {
    stop("ColData samples are not identical to colnames of expression data...",
         " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
  }
  pltData <- colDataNow
  
  # Calculate PCA
  if (! file.exists( "Data/bulkRNASeq/fullPCData_clusterOne.rda")) {
    pcDataFull <- prcomp(t(vsdHVG))
    pcData <- pcDataFull$x
    pcData <- pcData[,c(1:100)]
    save(pcData, file = "Data/bulkRNASeq/fullPCData_clusterOne.rda")
  } else {
    load("Data/bulkRNASeq/fullPCData_clusterOne.rda")
  }
  pltData$PC1_c1 <- pcData[,c(1)]
  pltData$PC2_c1 <- pcData[,c(2)]
  # Calculate Neighbors
  neighbors <- FindNeighbors(pcData)
  # Calculate clustering
  clusters <- FindClusters(neighbors$snn, resolution = .015)
  pltData$cluster_c1 <- clusters[,c(1)]
  # Calculate UMAP embedding  
  nnNow <- floor(length(colnames(vsdHVG)) * .2)
  umapData <- umap(t(vsdHVG), verbose = T,
                   n_neighbors = nnNow, 
                   min_dist = 1, pca = 100)
  pltData$UMAP_1_c1 <- umapData[,c(1)]
  pltData$UMAP_2_c1 <- umapData[,c(2)]
  save(pltData, file = "Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
} else {
  cat("\nUMAP data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
}
# Cluster 1 without the EWS genes
cat("\nRe-calculating VST for samples from cluster 1 with EWS geneset subtracted...\n")
if (! file.exists("Data/bulkRNASeq/nonEWSVSTCounts.rda")) {
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  load("Data/ewsGenes.rda")
  keepSamps <- clusterOneColData$samples
  nonEWSGenes <- rownames(expression)[! rownames(expression) %in% ewsGenes]
  expression <- expression[nonEWSGenes,keepSamps]
  # Make data Homoscedastic with VST
  vsd <- vst(expression, nsub = 1000)
  timestamp()
  cat("\nDone. Saving vst data...\n")
  save(vsd, file = "Data/bulkRNASeq/nonEWSVSTCounts_clusterOne.rda")
} else if (file.exists("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")) {
  cat("\nUMAP results found... Skipping this step...\n")
} else {
  cat("\nVST-transformed expression data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/nonEWSVSTCounts_clusterOne.rda")
}

cat("\nCalculating UMAP with EWS genes subtracted...\n")
if (! file.exists("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")) {
  # Calculate highly-variable genes (HVGs)
  rv <- matrixStats::rowVars(vsd)
  names(rv) <- rownames(vsd)
  rv <- rv[order(rv, decreasing = T)]
  hvgs <- rv[c(1:10000)]
  vsdHVG <- vsd[names(hvgs),]
  # dim(vsdHVG)
  rm(vsd)
  gc()
  colDataNow <- clusterOneColData[order(match(clusterOneColData$samples, colnames(vsdHVG))),]
  if (! all(colDataNow$samples == colnames(vsdHVG))) {
    stop("ColData samples are not identical to colnames of expression data...",
         " Please email Code/generateFigures_logFile.txt to author if you find this error and/or submit issue on github.")
  }
  pltData <- colDataNow
  
  # Calculate PCA
  if (! file.exists( "Data/bulkRNASeq/nonEWSPCData_clusterOne.rda")) {
    pcDataFull <- prcomp(t(vsdHVG))
    pcData <- pcDataFull$x
    pcData <- pcData[,c(1:100)]
    save(pcData, file = "Data/bulkRNASeq/nonEWSPCData_clusterOne.rda")
  } else {
    load("Data/bulkRNASeq/nonEWSPCData_clusterOne.rda")
  }
  pltData$PC1_c1 <- pcData[,c(1)]
  pltData$PC2_c1 <- pcData[,c(2)]
  # Calculate Neighbors
  neighbors <- FindNeighbors(pcData)
  # Calculate clustering
  clusters <- FindClusters(neighbors$snn, resolution = .015)
  pltData$cluster_c1 <- clusters[,c(1)]
  # Calculate UMAP embedding  
  nnNow <- floor(length(colnames(vsdHVG)) * .2)
  umapData <- umap(t(vsdHVG), verbose = T,
                   n_neighbors = nnNow, 
                   min_dist = 1, pca = 100)
  pltData$UMAP_1_c1 <- umapData[,c(1)]
  pltData$UMAP_2_c1 <- umapData[,c(2)]
  save(pltData, file = "Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")
} else {
  cat("\nUMAP data found! Loading and continuing to next step...\n")
  load("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")
}

cat("\nPlotting results of subclustering...\n")
load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
goodSamps <- pltData$samples[pltData$UMAP_2_c1 < 20] # Poor quality samples determined by UMAP
pltData <- pltData[pltData$samples %in% goodSamps,]
pltData$TissueFinal <- pltData$Tissue
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
pltData$Cluster <- pltData$cluster_c1
pltData$UMAP_1 <- pltData$UMAP_1_c1
pltData$UMAP_2 <- pltData$UMAP_2_c1
tableDF <- as.data.frame(table(pltData$Tissue, pltData$Cluster))
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Di.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Di_legend.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Dii.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Dii_legend.png", device = "png")
pltData$Tissue <- "Undefined"
pltData$Tissue[! pltData$TissueFinal %in% c("stem-like",
                                            "brain", "hESCs", 
                                            "Neural progenitor/stem cells",
                                            "fibroblasts", "iPSCs",  "MSCs",
                                            "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltData$Tissue[ pltData$TissueFinal %in% c("stem-like",
                                           "brain", "hESCs", 
                                           "Neural progenitor/stem cells",
                                           "fibroblasts", "iPSCs",  "MSCs",
                                           "ewing sarcoma", "neural crest cells")] <- pltData$TissueFinal[ pltData$TissueFinal %in% c("stem-like",
                                                                                                                                      "brain", "hESCs", 
                                                                                                                                      "Neural progenitor/stem cells",
                                                                                                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                      "ewing sarcoma", "neural crest cells")]
pltData$Tissue[pltData$Tissue %in% c("neural crest cells", "brain", "Neural progenitor/stem cells")] <- "Neural progenitors & brain"
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig1Diii.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Diii_legend.png", device = "png")
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resDF <- resDF[rowSums(resDF) > 100,]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(255)
colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")] <- stringr::str_to_sentence(colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")])
pheatmap(resMat, cluster_rows = F, cluster_cols = T, labels_col = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)
ap <- grid.grab()
plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig1Div.png", device = "png")
load("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")
pltData <- pltData[pltData$samples %in% goodSamps,]
pltData$TissueFinal <- pltData$Tissue
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
pltData$Cluster <- pltData$cluster_c1
pltData$UMAP_1 <- pltData$UMAP_1_c1
pltData$UMAP_2 <- pltData$UMAP_2_c1
tableDF <- as.data.frame(table(pltData$Tissue, pltData$Cluster))
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Dv.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Cluster")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 15,
       filename = "Figures/Fig1Dv_legend.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Dvi.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 15,
       filename = "Figures/Fig1Dvi_legend.png", device = "png")
pltData$Tissue <- "Undefined"
pltData$Tissue[! pltData$TissueFinal %in% c("stem-like",
                                            "brain", "hESCs", 
                                            "Neural progenitor/stem cells",
                                            "fibroblasts", "iPSCs",  "MSCs",
                                            "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltData$Tissue[ pltData$TissueFinal %in% c("stem-like",
                                           "brain", "hESCs", 
                                           "Neural progenitor/stem cells",
                                           "fibroblasts", "iPSCs",  "MSCs",
                                           "ewing sarcoma", "neural crest cells")] <- pltData$TissueFinal[ pltData$TissueFinal %in% c("stem-like",
                                                                                                                                      "brain", "hESCs", 
                                                                                                                                      "Neural progenitor/stem cells",
                                                                                                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                      "ewing sarcoma", "neural crest cells")]
pltData$Tissue[pltData$Tissue %in% c("neural crest cells", "brain", "Neural progenitor/stem cells")] <- "Neural progenitors & brain"
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend")
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig1Dvii.png", device = "png")
g1 <- ggplot(pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2",
                                           color = "Tissue")) +
  geom_point(size = .8) + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) 
ggsave(g1, height = 7.5, width = 15,
       filename = "Figures/Fig1Dvii_legend.png", device = "png")
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resDF <- resDF[rowSums(resDF) > 100,]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(255)
colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")] <- stringr::str_to_sentence(colnames(resMat)[! colnames(resMat) %in% c("HSCs", "iPSCs", "MSCs", "hESCs")])
pheatmap(resMat, cluster_rows = F, cluster_cols = T, labels_col = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)
ap <- grid.grab()
plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig1Dviii.png", device = "png")


cat("\nConstructing PHATE trajectories with and without EWS geneset...\n")
# Contstruct PHATE trajectories
kNow <- 560
# With EWS
if (! file.exists(paste0("Data/bulkRNASeq/PHATE_fullVST_nn", kNow, ".rda"))) {
  load("Data/bulkRNASeq/fullVSTCounts_clusterOne.rda")
  phateEWS <- phate(t(vsd), knn = kNow, n.jobs = 80, seed = 42, ndim = 3)
  phateEmbeddingEWS <- phateEWS$embedding
  save(phateEmbeddingEWS, file =  paste0("Data/bulkRNASeq/PHATE_fullVST_nn", kNow, ".rda"))
} else {
  cat("PHATE embedding found! Loading and continuing...")
  load(paste0("Data/bulkRNASeq/PHATE_fullVST_nn", kNow, ".rda"))
}
# Without EWS
if (! file.exists(paste0("Data/bulkRNASeq/PHATE_nonEWSVST_nn", kNow, ".rda"))) {
  load("Data/bulkRNASeq/nonEWSVSTCounts_clusterOne.rda")
  phateNonEWS <- phate(t(vsd), knn = kNow, n.jobs = 80, seed = 42, ndim = 3)
  phateEmbeddingNonEWS <- phateNonEWS$embedding
  save(phateEmbeddingNonEWS, file =  paste0("Data/bulkRNASeq/PHATE_nonEWSVST_nn", kNow, ".rda"))
} else {
  cat("PHATE embedding without EWS genes found! Loading and continuing...")
  load( paste0("Data/bulkRNASeq/PHATE_nonEWSVST_nn", kNow, ".rda"))
}

cat("\nCalculating DEGs and GSEA for PHATE trajectories in EWS cells & nonEWS-cells...")
if (! file.exists("Data/bulkRNASeq/trendyRes_nonEWSSamples.rda")) {
  load("Data/bulkRNASeq/fullRawCountsFiltered.rda")
  load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
  load(paste0("Data/bulkRNASeq/PHATE_fullVST_nn", kNow, ".rda"))
  pltData$PHATE_1 <- phateEmbeddingEWS[,c(1)]
  pltData$PHATE_2 <- phateEmbeddingEWS[,c(2)]
  pltData$PHATE_3 <- phateEmbeddingEWS[,c(3)]
  
  # All cells
  countsNow <- expression[,colnames(expression) %in% pltData$samples]
  Sizes <- MedianNorm(countsNow)
  countsNorm <- GetNormalizedMat(countsNow, Sizes)
  timeList <- list("PHATE_1" = range01(pltData$PHATE_1), 
                   "PHATE_2" = range01(pltData$PHATE_2), 
                   "PHATE_3" = range01(pltData$PHATE_3))
  resList <- lapply(timeList, FUN = function(timeNow) {
    doTrendyAnalysis(countsNorm, timeVec = timeNow)
  }) 
  save(resList, file = "Data/bulkRNASeq/trendyRes_allSamples.rda")
  resAllSamples <- resList
  
  # EWS cells
  ewsPhateColData <- pltData[pltData$Tissue == "ewing sarcoma",]
  countsNow <- expression[, colnames(expression) %in% ewsPhateColData$samples]
  Sizes <- MedianNorm(countsNow)
  countsNorm <- GetNormalizedMat(countsNow, Sizes)
  timeList <- list("PHATE_1" = range01(ewsPhateColData$PHATE_1), 
                   "PHATE_2" = range01(ewsPhateColData$PHATE_2), 
                   "PHATE_3" = range01(ewsPhateColData$PHATE_3))
  resList <- lapply(timeList, FUN = function(timeNow) {
    doTrendyAnalysis(countsNorm, timeVec = timeNow)
  }) 
  save(resList, file = "Data/bulkRNASeq/trendyRes_ewsSamples.rda")
  resEWSSamples <- resList
  
  # nonEWS cells
  nonEWSPhateColData <- pltData[pltData$Tissue != "ewing sarcoma",]
  countsNow <- expression[, colnames(expression) %in% nonEWSPhateColData$samples]
  Sizes <- MedianNorm(countsNow)
  countsNorm <- GetNormalizedMat(countsNow, Sizes)
  timeList <- list("PHATE_1" = range01(nonEWSPhateColData$PHATE_1), 
                   "PHATE_2" = range01(nonEWSPhateColData$PHATE_2), 
                   "PHATE_3" = range01(nonEWSPhateColData$PHATE_3))
  resList <- lapply(timeList, FUN = function(timeNow) {
    doTrendyAnalysis(countsNorm, timeVec = timeNow)
  }) 
  save(resList, file = "Data/bulkRNASeq/trendyRes_nonEWSSamples.rda")
  resNonEWSSamples <- resList
} else {
  cat("\nPhate DGE results found! Loading and continuing...\n")
  load("Data/bulkRNASeq/trendyRes_allSamples.rda")
  resAllSamples <- resList
  load("Data/bulkRNASeq/trendyRes_ewsSamples.rda")
  resEWSSamples <- resList
  load("Data/bulkRNASeq/trendyRes_nonEWSSamples.rda")
  resNonEWSSamples <- resList
}

cat("\nCalculating EWS-ETS fusion expression...\n")
if (! file.exists("../../../Preprocessing/RNA_Seq/EWS_Cells_EF1Quant/accList.txt")) {
  sraTab <- read.csv("Data/SRA_Accessions_SRR2GSM.csv", header = F, stringsAsFactors = F)
  ewsGSM <- pltData$samples[pltData$TissueType == "ewing sarcoma"]
  gsmSmall <- sraTab$V2
  gsmSmall <- gsub(gsmSmall, pattern = "(GSM[0-9]+)[\\.\\_]", replacement = "\\1")
  ewsInd <- grep(x = sraTab$V2, pattern = paste0(ewsGSM, collapse = "|"), perl = T)
  sraNow <- sraTab[ewsInd,]
  foundEWS <- gsub(sraNow$V2, pattern = "(GSM[0-9]+)_.+", replacement = "\\1")
  accList <- data.frame("Run" = sraNow$V1)
  write.table(accList, file = "../../../Preprocessing/RNA_Seq/EWS_Cells_EF1Quant/accList.txt", row.names = F, quote = F)
}

load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
pltData$PHATE_1 <- phateEmbeddingEWS[,c(1)]
pltData$PHATE_2 <- phateEmbeddingEWS[,c(2)]
pltData$PHATE_3 <- phateEmbeddingEWS[,c(3)]
pltData$Cluster <- pltData$cluster_c1
pltData$TissueFinal <- pltData$Tissue
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData, plotName = "Fig1Fi", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Cluster"))
DimPlot2(pltData, plotName = "Fig1Fii", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
pltData$Tissue <- "Undefined"
pltData$Tissue[! pltData$TissueFinal %in% c("stem-like",
                                            "brain", "hESCs", 
                                            "Neural progenitor/stem cells",
                                            "fibroblasts", "iPSCs",  "MSCs",
                                            "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltData$Tissue[ pltData$TissueFinal %in% c("stem-like",
                                           "brain", "hESCs", 
                                           "Neural progenitor/stem cells",
                                           "fibroblasts", "iPSCs",  "MSCs",
                                           "ewing sarcoma", "neural crest cells")] <- pltData$TissueFinal[ pltData$TissueFinal %in% c("stem-like",
                                                                                                                                      "brain", "hESCs", 
                                                                                                                                      "Neural progenitor/stem cells",
                                                                                                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                      "ewing sarcoma", "neural crest cells")]


pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData, plotName = "Fig1Fiii", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
pltData$specCat <- "Undefined"
pltData$specCat[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                           "Neural progenitor/stem cells")] <- pltData$Tissue[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                                                                                                         "Neural progenitor/stem cells")]
pltData2 <- pltData[pltData$specCat != "Undefined",]
pltData2$Tissue <- pltData2$specCat
DimPlot2(pltData2, plotName = "Fig1Fiv", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
p <- plot_ly(pltData, x = ~PHATE_1, size = .2, y = ~PHATE_2, 
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             z = ~PHATE_3, color = ~Tissue) %>%
  add_markers() 
Sys.setenv("plotly_username"="millerh1")
Sys.setenv("plotly_api_key"="s3xRqMxwBYfFnX3Au8d4")
api_create(p, filename = "Fig1Fv")
p <- plot_ly(pltData2, x = ~PHATE_1,
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             size = .2, y = ~PHATE_2, z = ~PHATE_3, 
             color = ~Tissue) %>%
  add_markers() 
api_create(p, filename = "Fig1Fvi")

load("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")
pltData$PHATE_1 <- phateEmbeddingNonEWS[,c(1)]
pltData$PHATE_2 <- phateEmbeddingNonEWS[,c(2)]
pltData$PHATE_3 <- phateEmbeddingNonEWS[,c(3)]
pltData$Cluster <- pltData$cluster_c1
pltData$TissueFinal <- pltData$Tissue
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData, plotName = "Fig1Fvii", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Cluster"))
DimPlot2(pltData, plotName = "Fig1Fviii", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
pltData$Tissue <- "Undefined"
pltData$Tissue[! pltData$TissueFinal %in% c("stem-like",
                                            "brain", "hESCs", 
                                            "Neural progenitor/stem cells",
                                            "fibroblasts", "iPSCs",  "MSCs",
                                            "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltData$Tissue[ pltData$TissueFinal %in% c("stem-like",
                                           "brain", "hESCs", 
                                           "Neural progenitor/stem cells",
                                           "fibroblasts", "iPSCs",  "MSCs",
                                           "ewing sarcoma", "neural crest cells")] <- pltData$TissueFinal[ pltData$TissueFinal %in% c("stem-like",
                                                                                                                                      "brain", "hESCs", 
                                                                                                                                      "Neural progenitor/stem cells",
                                                                                                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                      "ewing sarcoma", "neural crest cells")]


pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData, plotName = "Fig1Fix", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
pltData$specCat <- "Undefined"
pltData$specCat[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                           "Neural progenitor/stem cells")] <- pltData$Tissue[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                                                                                                         "Neural progenitor/stem cells")]
pltData2 <- pltData[pltData$specCat != "Undefined",]
pltData2$Tissue <- pltData2$specCat
DimPlot2(pltData2, plotName = "Fig1Fx", mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"))
p <- plot_ly(pltData, x = ~PHATE_1, size = .2, y = ~PHATE_2, 
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             z = ~PHATE_3, color = ~Tissue) %>%
  add_markers() 
api_create(p, filename = "Fig1Fxi")
p <- plot_ly(pltData2, x = ~PHATE_1,
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             size = .2, y = ~PHATE_2, z = ~PHATE_3, 
             color = ~Tissue) %>%
  add_markers() 
api_create(p, filename = "Fig1Fxii")

# Plot DGEs and GSEA from PHATE analysis
pltList <- lapply(names(resAllSamples), FUN = doPhateCorrPlot, resSamples = resAllSamples)
names(pltList) <- names(resAllSamples)
g1 <- ggarrange(plotlist = pltList, ncol = 3, nrow = 1)
ggsave(g1, height = 7.5, width = 24,
       filename = "Figures/Fig1Gi.png", device = "png")
pltList <- lapply(names(resEWSSamples), FUN = doPhateCorrPlot, resSamples = resEWSSamples)
names(pltList) <- names(resEWSSamples)
g1 <- ggarrange(plotlist = pltList, ncol = 3, nrow = 1)
ggsave(g1, height = 7.5, width = 24,
       filename = "Figures/Fig1Gii.png", device = "png")
pltList <- lapply(names(resNonEWSSamples), FUN = doPhateCorrPlot, resSamples = resNonEWSSamples)
names(pltList) <- names(resNonEWSSamples)
g1 <- ggarrange(plotlist = pltList, ncol = 3, nrow = 1)
ggsave(g1, height = 7.5, width = 24,
       filename = "Figures/Fig1Giii.png", device = "png")

# Do GSEA Plots
pltNames <- c("Miyagawa Targets of EWSR1-ETS Fusions Down" = "MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN",
              "Reactome Activation of ATR in Response to Replication Stress" = "REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS",
              "Hallmark Epithelial to Mesenchymal Transition" = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "PID Fanconi Pathway" = "PID_FANCONI_PATHWAY")
resSampsNow <- list(resAllSamples[[1]], resEWSSamples[[1]], resNonEWSSamples[[1]])
names(resSampsNow) <- c("All samples", "EWS samples", "Non-EWS samples")
pltList <- lapply(names(resSampsNow), FUN = function(nameNow, resSamples) {
  phateRes <- resSamples[[nameNow]]
  EGMTNow <- phateRes$EGMT
  pltListTmp <- lapply(names(pltNames), function(pathNow, pltNames, EGMT, nameNow) {
    pathPlot <- pltNames[[pathNow]]
    gseaPlot2(EGMT, title = pathNow,
              ID = pathPlot, scoreName = paste0(nameNow, " score"))
  }, pltNames = pltNames, EGMT = EGMTNow, nameNow = nameNow)
  ggarrange(plotlist = pltListTmp, ncol = 2, nrow = 2)
}, resSamples = resSampsNow)
g1 <- ggarrange(plotlist = pltList, ncol = 3, nrow = 1)
ggsave(g1, height = 15, width = 80, limitsize = F,
       filename = "Figures/Fig1Giv.png", device = "png")

# Compare EWS and Non-EWS
eresEWS <- resEWSSamples$PHATE_1$eres
eresNoEWS <- resNonEWSSamples$PHATE_1$eres
eresEWS <- eresEWS[order(abs(eresEWS$NES), decreasing = T),]
eresEWSTop500 <- eresEWS[c(1:500),]
eresNoEWS <- eresNoEWS[order(abs(eresNoEWS$NES), decreasing = T),]
eresNoEWSTop500 <- eresNoEWS[c(1:500),]
oList <- list(
  "Ewing sarcoma" = eresEWS$ID,
  "iPSCs/MSCs" = eresNoEWS$ID
)
plot.new()
dev.off()
gd <- venn.diagram(oList, filename = NULL, fill = c("firebrick", "skyblue"), 
                   margin = .05, cat.cex = 1.5, cat.dist = .05, cex = 1.4)
grid.draw(gd)
pt <- grid.grab()
gp1 <- as.ggplot(pt)
ggsave(gp1, height = 5, width = 7.3,
       filename = "Figures/Fig1Gv.png", device = "png")
eresEWSEx <- eresEWS[! eresEWS$ID %in% eresNoEWS$ID,]

# Exclusive EWS
GSEAResEWS <- resEWSSamples$PHATE_1
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "GO Cytosolic Ribosome",
                 ID = "GO_CYTOSOLIC_RIBOSOME", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxii.png", device = "png")
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "Reactome HATS Acetylate Histones",
                 ID = "REACTOME_HATS_ACETYLATE_HISTONES", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxiii.png", device = "png")
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "GO S-Adenosylmethionine-Dependent Methyltransferase Activity",
                 ID = "GO_S_ADENOSYLMETHIONINE_DEPENDENT_METHYLTRANSFERASE_ACTIVITY", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxiv.png", device = "png")
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "GO Error-free Translesional Synthesis",
                 ID = "GO_ERROR_FREE_TRANSLESION_SYNTHESIS", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxv.png", device = "png")
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "GO Mitochondrial Respiratory Chain Complex Assembly",
                 ID = "GO_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxvi.png", device = "png")
gp1 <- gseaPlot2(GSEAResEWS$EGMT, title = "GO Synapsis",
                 ID = "GO_SYNAPSIS", scoreName = "PHATE_1 score")
ggsave(gp1, height = 7.5, width = 12,
       filename = "Figures/Fig1Gxvii.png", device = "png")


# Generate plots for PHATE without EWS samples
kNow <- 560
# With EWS
if (! file.exists("Data/bulkRNASeq/PHATE_NoEWSCells_nn", kNow, ".rda")) {
  load("Data/bulkRNASeq/fullVSTCounts_clusterOne.rda")
  load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
  pltData <- pltData[pltData$Tissue != "ewing sarcoma",]
  vsd <- vsd[,colnames(vsd) %in% pltData$samples]
  all(pltData$samples == colnames(vsd))
  phateEWS <- phate(t(vsd), knn = kNow, n.jobs = 80, seed = 42, ndim = 3)
  phateEmbedding <- phateEWS$embedding
  save(phateEmbedding, file =  paste0("Data/bulkRNASeq/PHATE_NoEWSCells_nn", kNow, ".rda"))
} else {
  cat("Phate embedding found! Loading...")
  load(paste0("Data/bulkRNASeq/PHATE_NoEWSCells_nn", kNow, ".rda"))
}

load("Data/bulkRNASeq/fullUMAPData_clusterOne.rda")
pltData <- pltData[pltData$Tissue != "ewing sarcoma",]
pltData$PHATE_1 <- phateEmbedding[,c(1)]
pltData$PHATE_2 <- phateEmbedding[,c(2)]
pltData$PHATE_3 <- phateEmbedding[,c(3)]
pltData$Cluster <- pltData$cluster_c1
pltData$TissueFinal <- pltData$Tissue
pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Cluster"), plotName = "Fig1Hi", ptSize = .8)
DimPlot2(pltData2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"), plotName = "Fig1Hii", ptSize = .8)
pltData$Tissue <- "Undefined"
pltData$Tissue[! pltData$TissueFinal %in% c("stem-like",
                                            "brain", "hESCs", 
                                            "Neural progenitor/stem cells",
                                            "fibroblasts", "iPSCs",  "MSCs",
                                            "ewing sarcoma", "neural crest cells")] <- "Mixed/Other"
pltData$Tissue[ pltData$TissueFinal %in% c("stem-like",
                                           "brain", "hESCs", 
                                           "Neural progenitor/stem cells",
                                           "fibroblasts", "iPSCs",  "MSCs",
                                           "ewing sarcoma", "neural crest cells")] <- pltData$TissueFinal[ pltData$TissueFinal %in% c("stem-like",
                                                                                                                                      "brain", "hESCs", 
                                                                                                                                      "Neural progenitor/stem cells",
                                                                                                                                      "fibroblasts", "iPSCs",  "MSCs",
                                                                                                                                      "ewing sarcoma", "neural crest cells")]


pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")] <- stringr::str_to_sentence(pltData$Tissue[! pltData$Tissue %in% c("hESCs", "HSCs", "MSCs", "iPSCs")])
DimPlot2(pltData, mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"), plotName = "Fig1Hiii", ptSize = .8)
pltData$specCat <- "Undefined"
pltData$specCat[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                           "Neural progenitor/stem cells")] <- pltData$Tissue[pltData$TissueFinal %in% c("ewing sarcoma", "HSCs", "MSCs", "iPSCs", "neural crest cells",
                                                                                                                         "Neural progenitor/stem cells")]
pltData2 <- pltData[pltData$specCat != "Undefined",]
pltData2$Tissue <- pltData2$specCat
DimPlot2(pltData2, mapping = aes_string(x = "PHATE_1", y = "PHATE_2", color = "Tissue"), plotName = "Fig1Hiv", ptSize = .8)
p <- plot_ly(pltData, x = ~PHATE_1, size = .2, y = ~PHATE_2, 
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             z = ~PHATE_3, color = ~Tissue) %>%
  add_markers() 
Sys.setenv("plotly_username"="millerh1")
Sys.setenv("plotly_api_key"="s3xRqMxwBYfFnX3Au8d4")
api_create(p, filename = "Fig1Hv")
p <- plot_ly(pltData2, x = ~PHATE_1,
             colors = colorRampPalette(brewer.pal(10,"Spectral"))(255),
             size = .2, y = ~PHATE_2, z = ~PHATE_3, 
             color = ~Tissue) %>%
  add_markers() 
api_create(p, filename = "Fig1Hvi")

# Comparison plots
pltData$PHATE_EWS <- pltData$PHATE_1
phateDataNow <- pltData
phateDataNow <- phateDataNow[,c(1, 22, 24)]

load( paste0("Data/bulkRNASeq/PHATE_nonEWSVST_nn", kNow, ".rda"))
load("Data/bulkRNASeq/nonEWSUMAPData_clusterOne.rda")
pltData$PHATE_1 <- phateEmbeddingNonEWS[,c(1)]
pltData$PHATE_2 <- phateEmbeddingNonEWS[,c(2)]
pltData$PHATE_3 <- phateEmbeddingNonEWS[,c(3)]
pltData$Cluster <- pltData$cluster_c1
pltData$TissueFinal <- pltData$Tissue
pltData$PHATE_Normal <- pltData$PHATE_1
phateDataNow2 <- pltData[,c(1, 22, 23)]
phateDataNowFinal <- merge(x = phateDataNow, y =phateDataNow2,
                           by = c("samples", "TissueFinal"))
phateDataNowFinal <- phateDataNowFinal[phateDataNowFinal$TissueFinal %in% c("iPSCs", "MSCs"),]
phateDataNowFinal$Tissue <- phateDataNowFinal$TissueFinal
DimPlot2(phateDataNowFinal, mapping = aes_string(x = "PHATE_EWS", y = "PHATE_Normal", color = "Tissue"), plotName = "Fig1Hvii")

#####################################################################################
################################ Single Cell RNASeq #################################
#####################################################################################

# Part I of single cell -- PHATE_1 signature in EWS datasets
if (! file.exists("Data/scRNASeq/BDInt.rda")) {
  # Add the Bishop group's single cell
  # -- CHLA9
  CHLA9 <- Read10X("Data/scRNASeq/nonPDB_matrices/CHLA9/")
  genesFinal <- rownames(CHLA9)
  CHLA9 <- CreateSeuratObject(CHLA9)
  CHLA9$Sample <- "CHLA9"
  # CHLA10
  CHLA10 <- Read10X("Data/scRNASeq/nonPDB_matrices/CHLA10/")
  genesFinal <- rownames(CHLA10)
  CHLA10 <- CreateSeuratObject(CHLA10)
  CHLA10$Sample <- "CHLA10"
  # TC71
  TC71 <- Read10X("Data/scRNASeq/nonPDB_matrices/TC71/")
  genesFinal <- rownames(TC71)
  TC71 <- CreateSeuratObject(TC71)
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
         filename = "Figures/Fig2Ai.png", device = "png")
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
         filename = "Figures/Fig2Aii.png", device = "png")
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
         filename = "Figures/Fig2Aiii.png", device = "png")
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
         filename = "Figures/Fig2Aiv.png", device = "png")
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
         filename = "Figures/Fig2Av.png", device = "png")
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
         filename = "Figures/Fig2Avi.png", device = "png")
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
  
  # Add in the Delattre group's single cell
  if (! file.exists("Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX/GSM3730317_PDX-861_matrix.mtx.gz")) {
    download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130025/suppl/GSE130025_RAW.tar", destfile = "Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX.tar")
    untar("Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX.tar", exdir = "Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX")
    rmFiles <- list.files("Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX/", pattern = "\\.txt|\\.bed", full.names = T)
    file.remove(rmFiles)
  }
  filesNow <- list.files("Data/scRNASeq/nonPDB_matrices/Delattre_EWS_PDX", pattern = "mtx", full.names = T)
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
         filename = "Figures/Fig2Bi.png", device = "png")
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
         filename = "Figures/Fig2Bii.png", device = "png")
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
         filename = "Figures/Fig2Biii.png", device = "png")
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
         filename = "Figures/Fig2Biv.png", device = "png")
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
         filename = "Figures/Fig2Bv.png", device = "png")
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
         filename = "Figures/Fig2Bvi.png", device = "png")
  
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
  BDInt <- IntegrateData(anchorset = srtList, dims = 1:30)
  DefaultAssay(BDInt) <- "integrated"
  BDInt <- ScaleData(BDInt)
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
  pltData$Outlier[pltData$UMAP_1 < -4.5 & pltData$UMAP_2 > .3] <- T
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
  BDInt <- RunUMAP(BDInt, dims = 1:30)
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

# Cluster heatmap
tableDF <- as.data.frame(table(pltData$SRS, pltData$Cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
pheatmap(t(resMat), cluster_rows = T, cluster_cols = F, labels_row = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 18)
ap <- grid.grab()
plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig2Cxi.png", device = "png")
tableDF <- as.data.frame(table(pltData$Phase, pltData$Cluster))
resDF <- tableDF %>% spread(key = Var1, value = Freq)
rownames(resDF) <- paste0("Cluster_", resDF$Var2)
resDF <- resDF[,c(-1)]
resMat <- apply(resDF, MARGIN = 2, FUN = scale, center = F, scale = T)
rownames(resMat) <- rownames(resDF)
paletteLength <- 255
heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
pheatmap(t(resMat), cluster_rows = T, cluster_cols = F, labels_row = as.character(colnames(resMat)),
         color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 18)
ap <- grid.grab()
plt <- ggplotify::as.ggplot(ap)
ggsave(plt, height = 5.5, width = 14,
       filename = "Figures/Fig2Cxii.png", device = "png")

if (! file.exists("Data/scRNASeq/BDIntMarkers.rda")) {
  # Find cluster markers
  BDIntMarkers <- FindAllMarkers(BDInt, logfc.threshold = .125)
  save(BDIntMarkers,  file = "Data/scRNASeq/BDIntMarkers.rda")
  doClusterMarkerPlot(srt = BDInt, srtMarkers = BDIntMarkers, outName = "Fig2Cxiii")
} else {
  load("Data/scRNASeq/BDIntMarkers.rda")  
}

# Plot PHATE_1 effect from bulk RNA-Seq
DefaultAssay(BDInt) <- "RNA"
BDInt <- NormalizeData(BDInt)  
BDInt <- ScaleData(BDInt, do.center = F, 
                   block.size = 10000)
load("Data/bulkRNASeq/trendyRes_ewsSamples.rda")
rankDF <- resEWSSamples$PHATE_1$rankDF 
nNow <- 100
topRankUp <- rankDF %>% top_n(n = nNow, wt = value)
topPhateOneUp <- topRankUp$geneName
countsNow <- BDInt@assays$RNA@scale.data[rownames(BDInt@assays$RNA@scale.data)
                                         %in% topPhateOneUp,]
phate_1_score <- colMeans(countsNow)
BDInt$PHATE_1_UP_SCORE <- phate_1_score
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
g1 <- ggplot(pltData, mapping = aes_string(x = "PHATE_1_UP_SCORE",
                                           y = "Cluster",
                                           fill = "Cluster")) +
  geom_density_ridges(scale = 2) +
  theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  rremove("legend") + rremove("ylab") + 
  xlab("PHATE_1 score")
ggsave(g1, filename = "Figures/Fig2Cxiv.png", height = 7.5, width = 9)


# Part II of single cell -- integration of normal and EWS datasets
if (! file.exists("Data/scRNASeq/NormalCellsInt.rda")) {
  # Integrate with iPSCs, MSCs, NCCs, NPCs, etc
  keepStudies <- c("SRA710104", "SRA795539", "SRA843432")
  srtFinalList <- list()
  for (i in 1:length(keepStudies)) {
    studyNow <- keepStudies[i]
    print(studyNow)
    sms <- list.files("Data/scRNASeq/var/www/html/SRA/SRA.final/",
                      pattern = studyNow, full.names = T)
    names(sms) <- gsub(sms, pattern = ".+SRA[0-9]+_(SRS[0-9]+).sparse.RData", 
                       replacement = "\\1")
    srtList <- list()
    for (j in 1:length(sms)) {
      SRSNow <- names(sms)[j]
      print(SRSNow)
      load(sms[j])
      # -- standardize row names
      genes <- rownames(sm)
      genesNow <- gsub(genes, pattern = "(.+)_ENSG.+", replacement = "\\1")
      badInds <- which(duplicated(genesNow))
      genesFinal <- genesNow[-badInds]
      smFinal <- sm[-badInds,]
      rownames(smFinal) <- genesFinal
      # -- make seurat object
      srtNow <- CreateSeuratObject(smFinal)
      srtNow$SRA <- studyNow
      srtNow$SRS <- SRSNow
      srtList[[j]] <- srtNow
    }
    srtTmp <- merge(x = srtList[[1]], y = srtList[c(2:length(srtList))])
    srtFinalList[[i]] <- srtTmp
  }
  names(srtFinalList) <- keepStudies
  
  ## QC each one
  # SRA710104
  srtFinalList[["SRA710104"]][["percent.mt"]] <- PercentageFeatureSet(srtFinalList[["SRA710104"]], 
                                                                      pattern = "^MT-")
  pltData <- srtFinalList[["SRA710104"]]@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
    rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Di.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Dii.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Diii.png", device = "png")
  # QC each
  srtFinalList[["SRA710104"]] <- autoQCSRT(srtFinalList[["SRA710104"]], minNumCounts = 3500,
                                           whichColumn = "SRS",
                                           minNumGenes = 600, maxPercentMT = 20 )
  # Show again after QC
  pltData <- srtFinalList[["SRA710104"]]@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Div.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Dv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 12,
         filename = "Figures/Fig2Dvi.png", device = "png")
  
  # SRA795539
  srtFinalList[["SRA795539"]][["percent.mt"]] <- PercentageFeatureSet(srtFinalList[["SRA795539"]], 
                                                                      pattern = "^MT-")
  pltData <- srtFinalList[["SRA795539"]]@meta.data
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
         filename = "Figures/Fig2Ei.png", device = "png")
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
         filename = "Figures/Fig2Eii.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 8,
         filename = "Figures/Fig2Eiii.png", device = "png")
  # QC each
  srtFinalList[["SRA795539"]] <- autoQCSRT(srtFinalList[["SRA795539"]], 
                                           minNumCounts = 7500,
                                           whichColumn = "SRS",
                                           minNumGenes = 2000, maxPercentMT = 15 )
  # Show again after QC
  pltData <- srtFinalList[["SRA795539"]]@meta.data
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
         filename = "Figures/Fig2Eiv.png", device = "png")
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
         filename = "Figures/Fig2Ev.png", device = "png")
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
         filename = "Figures/Fig2Evi.png", device = "png")
  # SRA843432
  srtFinalList[["SRA843432"]][["percent.mt"]] <- PercentageFeatureSet(srtFinalList[["SRA843432"]], 
                                                                      pattern = "^MT-")
  pltData <- srtFinalList[["SRA843432"]]@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fi.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fii.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .025) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fiii.png", device = "png")
  # QC each
  srtFinalList[["SRA843432"]] <- autoQCSRT(srtFinalList[["SRA843432"]], minNumCounts = 1500,
                                           whichColumn = "SRS",
                                           minNumGenes = 400, maxPercentMT = 15 )
  # Show again after QC
  pltData <- srtFinalList[["SRA843432"]]@meta.data
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "percent.mt",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Mitochrondrial counts (%)") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fiv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nFeature_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of genes expressed") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fv.png", device = "png")
  g1 <- ggplot(pltData, mapping = aes_string(x = "SRS", y = "nCount_RNA",
                                             color = "SRS")) +
    geom_violin(size = .8) + theme_pubr(border = T, base_size = 22) +
    guides(colour = guide_legend(override.aes = list(size=3))) + 
    theme(legend.position="right") +
    ylab("Number of read counts") + geom_jitter(alpha = .05) +
    theme(#border = element_line(colour = "black", size = 1),
      panel.border = element_rect(colour = "black", fill=NA, size=1)) + rremove("legend") +
    rremove("xlab")
  ggsave(g1, height = 7.5, width = 16,
         filename = "Figures/Fig2Fvi.png", device = "png")
  
  # Merge together
  srtFinal <- merge(x = srtFinalList[[1]], y = srtFinalList[c(2:length(srtFinalList))])
  srtFinal$Sample <- rownames(srtFinal@meta.data)
  studyTissueTable <- data.frame(
    SRS = unique(srtFinal$SRS),
    Tissue = c("aMSCs", "aMSCs", "aMSCs", 
               "hNCCs", "hNCCs", 
               "iPSCs", "iPSCs", "iPSCs", "NPCs", "NPCs", "NPCs")
  )
  srtMeta <- srtFinal@meta.data
  srtMetaMerge <- merge(x = srtMeta, y = studyTissueTable, by = "SRS")
  srtMetaMerge <- srtMetaMerge[order(match(srtMetaMerge$Sample, srtFinal$Sample)),]
  srtFinal$Tissue <- srtMetaMerge$Tissue
  srtFinal <- downsampleSRT(srtFinal, n = 4000, sampleCol = "Tissue")
  normalCells <- srtFinal
  
  # Normal clustering
  normalCells <- NormalizeData(normalCells)
  normalCells <- FindVariableFeatures(normalCells)
  normalCells <- ScaleData(normalCells)
  normalCells <- RunPCA(normalCells)
  normalCells <- FindNeighbors(normalCells)
  normalCells <- FindClusters(normalCells)
  normalCells <- RunUMAP(normalCells, dims = 1:50)
  normalCells <- CellCycleScoring(normalCells, 
                                  s.features = cc.genes$s.genes,
                                  g2m.features = cc.genes$g2m.genes)
  pltData <- normalCells@meta.data
  pltData$UMAP_1 <- normalCells@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- normalCells@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Sample"), 
           plotName = "Fig2Gi")
  # Needs batch correct 
  srtList <- SplitObject(normalCells, split.by = "SRS")
  for (i in 1:length(srtList)) {
    srtList[[i]] <- NormalizeData(srtList[[i]], verbose = FALSE)
    srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                         nfeatures = 5000, verbose = FALSE)
  }
  srtList <- FindIntegrationAnchors(object.list = srtList, dims = 1:30,
                                    anchor.features = 5000)
  normalCellsInt <- IntegrateData(anchorset = srtList, dims = 1:30)
  DefaultAssay(normalCellsInt) <- "integrated"
  normalCellsInt <- ScaleData(normalCellsInt)
  normalCellsInt <- RunPCA(normalCellsInt)
  normalCellsInt <- FindNeighbors(normalCellsInt)
  normalCellsInt <- FindClusters(normalCellsInt)
  normalCellsInt <- RunUMAP(normalCellsInt, dims = 1:30)
  
  # Separate plots
  pltData <- normalCellsInt@meta.data
  pltData$UMAP_1 <- normalCellsInt@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- normalCellsInt@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  pltData$Cluster <- pltData$seurat_clusters
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"), 
           plotName = "Fig2Gii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Tissue"), 
           plotName = "Fig2Giii")
  DimPlot2(pltData, 
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"), 
           plotName = "Fig2Giv")
  
  save(normalCellsInt, file = "Data/scRNASeq/NormalCellsInt.rda")
} else {
  cat("\nIntegrated normal tissues object found! Loading and continuing...")
  load("Data/scRNASeq/NormalCellsInt.rda")
  pltData <- normalCellsInt@meta.data
  pltData$UMAP_1 <- normalCellsInt@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- normalCellsInt@reductions$umap@cell.embeddings[,c(2)]
  pltData$Sample <- pltData$SRS
  pltData$Cluster <- pltData$seurat_clusters
}


# Generate figures for single cell CCA in EWS, iPSCs, and MSCs
if (! file.exists("Data/scRNASeq/srtWithEWS_iPSC_MSC_EWS_Integrated.rda")) {
  # Compare EWS with Cluster One Cell types
  if (! file.exists("Data/scRNASeq/srtIntEWS.rda")) {
    BDInt$Tissue <- BDInt$SRS
    BDInt$Tissue[grep(BDInt$SRS, pattern = "PDX")] <- "EWS PDX"
    BDInt$Tissue[BDInt$SRA == "BISHOP"] <- "EWS Cells"
    srt <- merge(x = normalCellsInt, y = BDInt)
    keep <- c("iPSCs", "aMSCs", "DELATTRE", "BISHOP")
    keepCells <- rownames(srt@meta.data)[srt@meta.data$Tissue %in% keep |
                                           srt@meta.data$SRA %in% keep]
    srtNow <- subset(srt, cells = keepCells)
    srtNow <- downsampleSRT(srtNow, n = 2000, sampleCol = "Tissue")
    srtList <- SplitObject(srtNow, split.by = "SRS")
    ccGenes <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))
    for (i in 1:length(srtList)) {
      print(i)
      DefaultAssay(srtList[[i]]) <- "RNA"
      srtList[[i]] <- NormalizeData(srtList[[i]], verbose = T)
      srtList[[i]] <- CellCycleScoring(srtList[[i]], s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
      srtList[[i]] <- ScaleData(object = srtList[[i]], block.size = 10000, 
                                vars.to.regress = c("S.Score", "G2M.Score"))
      srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                           nfeatures = 5000, verbose = T)
      srtList[[i]]@assays$RNA@var.features <- srtList[[i]]@assays$RNA@var.features[! srtList[[i]]@assays$RNA@var.features %in% ccGenes]
    }
    srtAnchors <- FindIntegrationAnchors(object.list = srtList, 
                                         scale = FALSE, anchor.features = 4500,
                                         verbose = T, dims = 1:30)
    srtInt <- IntegrateData(anchorset = srtAnchors, dims = 1:30)
    DefaultAssay(srtInt) <- "integrated"
    srtInt <- ScaleData(srtInt, block.size = 10000)
    srtInt <- RunPCA(srtInt, npcs = 50)
    srtIntEWS <- srtInt
    srtIntEWS <- FindNeighbors(srtIntEWS)
    srtIntEWS <- FindClusters(srtIntEWS)
    srtIntEWS <- RunUMAP(srtIntEWS, dims = 1:50)
    save(srtIntEWS, file = "Data/scRNASeq/srtIntEWS.rda")
  } else {
    cat("\nIntegrated EWS & Normal cells object found! Loading and continuing...\n")
    load("Data/scRNASeq/srtIntEWS.rda")
  }
  
  # Without EWS Genes
  if (! file.exists("Data/scRNASeq/srtIntNonEWS.rda")) {
    load("Data/ewsGenes.rda")
    for (i in 1:length(srtList)) {
      print(i)
      DefaultAssay(srtList[[i]]) <- "RNA"
      srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                           nfeatures = 6000, verbose = T)
      srtList[[i]]@assays$RNA@var.features <- srtList[[i]]@assays$RNA@var.features[! srtList[[i]]@assays$RNA@var.features %in% ccGenes &
                                                                                     ! srtList[[i]]@assays$RNA@var.features %in% ewsGenes]
      print(length(srtList[[i]]@assays$RNA@var.features))
    }
    srtAnchors <- FindIntegrationAnchors(object.list = srtList, 
                                         scale = FALSE, anchor.features = 4500,
                                         verbose = T, dims = 1:30)
    srtInt <- IntegrateData(anchorset = srtAnchors, dims = 1:30)
    DefaultAssay(srtInt) <- "integrated"
    srtInt <- ScaleData(srtInt, block.size = 10000)
    srtInt <- RunPCA(srtInt, npcs = 50)
    srtIntNonEWS <- srtInt
    srtIntNonEWS <- FindNeighbors(srtIntNonEWS)
    srtIntNonEWS <- FindClusters(srtIntNonEWS)
    srtIntNonEWS <- RunUMAP(srtIntNonEWS, dims = 1:50)
    save(srtIntNonEWS, file = "Data/scRNASeq/srtIntNonEWS.rda") 
  } else {
    cat("\nIntegrated EWS & Normal cells object without EWS genes found! Loading and continuing...\n")
    load("Data/scRNASeq/srtIntNonEWS.rda")
  }
   
  # Plot results with EWS genes
  pltData <- srtIntEWS@meta.data
  pltData$UMAP_1 <- srtIntEWS@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- srtIntEWS@reductions$umap@cell.embeddings[,c(2)]
  pltData$Cluster <- pltData$seurat_clusters
  DimPlot2(pltData, plotName = "Fig2Hi",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Tissue"))
  DimPlot2(pltData, plotName = "Fig2Hii",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"))
  DimPlot2(pltData, plotName = "Fig2Hiii",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"))
  tableDF <- as.data.frame(table(pltData$Tissue, pltData$seurat_clusters))
  resDF <- tableDF %>% spread(key = Var1, value = Freq)
  rownames(resDF) <- paste0("Cluster_", resDF$Var2)
  resDF <- resDF[,c(-1)]
  resMat <- apply(t(resDF), MARGIN = 2, FUN = scale, center = F, scale = T)
  rownames(resMat) <- colnames(resDF)
  paletteLength <- 255
  heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
  pheatmap(t(resMat), cluster_rows = F, cluster_cols = T, labels_row = as.character(colnames(resMat)),
           color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)
  ap <- grid.grab()
  plt <- ggplotify::as.ggplot(ap)
  ggsave(plt, height = 5.5, width = 14,
         filename = "Figures/Fig2Hiv.png", device = "png")
  
  # Non-EWS
  pltData <- srtIntNonEWS@meta.data
  pltData$UMAP_1 <- srtIntNonEWS@reductions$umap@cell.embeddings[,c(1)]
  pltData$UMAP_2 <- srtIntNonEWS@reductions$umap@cell.embeddings[,c(2)]
  pltData$Cluster <- pltData$seurat_clusters
  DimPlot2(pltData, plotName = "Fig2Hv",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Tissue"))
  DimPlot2(pltData, plotName = "Fig2Hvi",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"))
  DimPlot2(pltData, plotName = "Fig2Hvii",
           mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"))
  tableDF <- as.data.frame(table(pltData$Tissue, pltData$seurat_clusters))
  resDF <- tableDF %>% spread(key = Var1, value = Freq)
  rownames(resDF) <- paste0("Cluster_", resDF$Var2)
  resDF <- resDF[,c(-1)]
  resMat <- apply(t(resDF), MARGIN = 2, FUN = scale, center = F, scale = T)
  rownames(resMat) <- colnames(resDF)
  paletteLength <- 255
  heatColor <- colorRampPalette((brewer.pal(9, "Reds")))(paletteLength)
  pheatmap(t(resMat), cluster_rows = F, cluster_cols = T, labels_row = as.character(colnames(resMat)),
           color = heatColor, angle_col = 315, margins = c(5, 19), fontsize = 13.5)
  ap <- grid.grab()
  plt <- ggplotify::as.ggplot(ap)
  ggsave(plt, height = 5.5, width = 14,
         filename = "Figures/Fig2Hviii.png", device = "png")
  
} else {
  cat("Integrated SRT with EWS genes found! Loading...")
  load("Data/scRNASeq/srtWithEWS_iPSC_MSC_EWS_Integrated.rda")
}


if (! file.exists("Data/scRNASeq/srtWithOutEWS_iPSC_MSC_EWS_Integrated.rda")) {
  # Without EWS genes
  load("Data/scRNASeq/srtHuman_nonEWS_Processed.rda")
  srt <- srtNonEWS
  keep <- c("MSCs", "iPSCs",  "Ewing sarcoma")
  keepCells <- rownames(srt@meta.data)[srt@meta.data$specificTissueType %in% keep]
  srtNow <- subset(srt, cells = keepCells)
  srtList <- SplitObject(srtNow, split.by = "SRS")
  for (i in 1:length(srtList)) {
    srtList[[i]] <- NormalizeData(srtList[[i]], verbose = T)
    srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = T)
  }
  srtAnchors <- FindIntegrationAnchors(object.list = srtList, verbose = T, dims = 1:30)
  srtInt <- IntegrateData(anchorset = srtAnchors, dims = 1:30)
  DefaultAssay(srtInt) <- "integrated"
  srtInt <- ScaleData(srtInt)
  srtInt <- RunPCA(srtInt, npcs = 50)
  srtIntNonEWS <- srtInt
  srtIntNonEWS <- RunUMAP(srtIntNonEWS, dims = 1:40, min.dist = .05)
  save(srtIntNonEWS, file = "Data/scRNASeq/srtWithOutEWS_iPSC_MSC_EWS_Integrated.rda")
} else {
  cat("Integrated SRT without EWS genes found! Loading...")
  load("Data/scRNASeq/srtWithOutEWS_iPSC_MSC_EWS_Integrated.rda")
}

pltData <- srtIntEWS@meta.data
pltData$UMAP_1 <- srtIntEWS@reductions$umap@cell.embeddings[,c(1)]
pltData$UMAP_2 <- srtIntEWS@reductions$umap@cell.embeddings[,c(2)]
pltData$Tissue <- pltData$specificTissueType
DimPlot2(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Tissue"), plotName = "Fig2Hi")
DimPlot2(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "SRS"), plotName = "Fig2Hii")
DimPlot2(data = pltData, mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"), plotName = "Fig2Hiii")


dp1 <- DimPlot(srtIntEWS, reduction = "umap",pt.size = .6,
               group.by = "specificTissueType")
g1 <- dp1 + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  rremove("legend") 
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig2Hi.png", device = "png")
g1 <- dp1 + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") + labs(color = "Tissue") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig2Ci_legend.png", device = "png")

dp1 <- DimPlot(srtIntNonEWS, reduction = "umap",pt.size = .6,
               group.by = "specificTissueType")
g1 <- dp1 + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  rremove("legend") 
ggsave(g1, height = 7.5, width = 8,
       filename = "Figures/Fig2Cii.png", device = "png")
g1 <- dp1 + theme_pubr(border = T, base_size = 22) +
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="right") + labs(color = "Tissue") +
  theme(#border = element_line(colour = "black", size = 1),
    panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave(g1, height = 7.5, width = 12,
       filename = "Figures/Fig2Cii_legend.png", device = "png")


# Imput bulk RNA-Seq PHATE_1 into single cell






# Merge EWS and Normal Tissues

BDIntNow <- downsampleSRT(BDInt, n = 500, sampleCol = "SRS")
DefaultAssay(BDIntNow) <- "RNA"
DefaultAssay(normalCellsInt) <- "RNA"
goodCells <- rownames(normalCellsInt@meta.data)[normalCellsInt$Tissue %in% c("aMSCs", "iPSCs")]
normalCellsIntNow <- subset(normalCellsInt, cells = goodCells)
normalCellsIntNow <- downsampleSRT(normalCellsIntNow, n = 700, sampleCol = "SRS")
mergeSrtNow <- merge(x = BDIntNow, y = normalCellsIntNow)
# Needs batch correct 
srtList <- SplitObject(mergeSrtNow, split.by = "SRS")
for (i in 1:length(srtList)) {
  print(i)
  srtList[[i]] <- NormalizeData(srtList[[i]], verbose = FALSE)
  srtList[[i]] <- FindVariableFeatures(srtList[[i]], selection.method = "vst", 
                                       nfeatures = 5000, verbose = FALSE)
}
srtList <- FindIntegrationAnchors(object.list = srtList, dims = 1:30,
                                  #eps = .5,
                                  anchor.features = 5000)
fullInt <- IntegrateData(anchorset = srtList, dims = 1:30)
DefaultAssay(fullInt) <- "integrated"
fullInt <- ScaleData(fullInt)
fullInt <- RunPCA(fullInt)
fullInt <- FindNeighbors(fullInt)
fullInt <- FindClusters(fullInt)
fullInt <- RunUMAP(fullInt, dims = 1:30)
fullInt$Tissue[fullInt$SRS %in% c("CHLA10", "CHLA9", "TC71")] <- "EWS"
fullInt$Tissue[fullInt$SRS %in% c("PDX-1058", "PDX-856","PDX-861",
                                  "PDX-184", "PDX-352")] <- "EWS"

# Separate plots
pltData <- fullInt@meta.data
pltData$UMAP_1 <- fullInt@reductions$umap@cell.embeddings[,c(1)]
pltData$UMAP_2 <- fullInt@reductions$umap@cell.embeddings[,c(2)]
pltData$Sample <- pltData$SRS
pltData$Cluster <- pltData$seurat_clusters
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Cluster"), 
         plotName = "Fig2Hii")
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Tissue"), 
         plotName = "Fig2Hiii")
DimPlot2(pltData, 
         mapping = aes_string(x = "UMAP_1", y = "UMAP_2", color = "Phase"), 
         plotName = "Fig2Hiv")

save(fullInt, file = "Data/scRNASeq/fullInt.rda")



