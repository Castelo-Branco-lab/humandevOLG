# ArchR pipeline for snATAC analysis 
# Load libraries
library(ArchR)
library(Seurat)
library(Signac)
library(plyr)
library(ggplot2)
library(cowplot)
library(hues)
set.seed(1)
addArchRThreads(threads = 1)
addArchRGenome("hg38")

# Get CellRanger ATAC output and create arrow files
sample_ids <- c("HES1", "P18861_1001", "P18861_1002")
input_files <- paste0("/data/proj/GCB_FBP/human_dev/data/scatac/cellranger/", sample_ids, "/outs/fragments.tsv.gz")
names(input_files) <- sample_ids

# Create arrow files - filter cells with min TSS enrichment > 10 and min number of Fragments > 1000
# Genomic windows - 2kb
ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  minTSS = 10, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list("tileSize" = 2000)
)

# Group samples in ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "hca_dev_best3",
  copyArrows = TRUE
)

# Dimensionality reduction with IterativeLSI
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 3, 
    clusterParams = list(
        resolution = 0.3,
        maxClusters = 8,
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    corCutOff = 0.8,
    depthCol = "nFrags"
)    

# Clutering
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", force = TRUE, resolution = 1)

# UMAP representation
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)

# UMAP plots before integration
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-Tile.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Integration with Hamorny
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

# UMAP plots after integration
p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
plotPDF(p, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Integrate with RNA data
load("/data/proj/GCB_DvB/humanPublication/SummarizedDatasethFetalrna.Rdata")
ematu <- networkExpressionFile[grepl("_u",row.names(networkExpressionFile)),]
rownames(ematu) <- gsub("_u", "", rownames(ematu))
ann <- networkAnnotableFile[colnames(ematu), ]

seurat_rna <- CreateSeuratObject(counts=ematu)
AddMetaData(seurat_rna, networkAnnotableFile) 
seurat_rna$seurat_clusters <- networkAnnotableFile$seurat_clusters
seurat_rna$version <- networkAnnotableFile$version

## Add cluster names
new.cluster.ids <- c("Excitatory neurons cortex",
                     "Radial Glia potential glioblast",
                     "Inhibitory neurons Midbrain, possibly GABAergic",
                     "Excitatory neurons midbrain",
                     "Forebrain early neuroblast possibly GABAergic", 
                     "Cortical Interneurons",
                     "Excitatory neurons possibly midbrain", 
                     "Glioblast", 
                     "Radial Glia/Glioblast", 
                     "Cortical Pyramidal", 
                     "Forebrain inhibitory neuroblast",
                     "Radial Glia cycling",
                     "GABAergic forebrain",
                     "Mid/Hindbrain neuroblast",
                     "Forebrain neural progenitor EMX1",
                     "Glioblast/Pre-OPC",
                     "Neuroblast motorneuron/GABAergic?",
                     "Forebrain neural progenitor EMX1",
                     "Radial Glia",
                     "U", 
                     "Radial Glia VLMC primed?", 
                     "Hindbrain neuroblast",
                     "GABAergic or interneuron neuroblast probably midbrain",
                     "Radial Glia",
                     "Striatum/Cortical neurons",
                     "Radial glia/Glioblast/Forebrain progenitor EMX1",
                     "OPCs",
                     "VLMCs",
                     "Midbrain inhibitory neuroblast", 
                     "Endothelial",
                     "Microglia",
                     "Midbrain/Hindbrain inhibitory neuroblast")

names(new.cluster.ids) <- levels(seurat_rna$seurat_clusters)
seurat_rna$annotated_clusters <- mapvalues(seurat_rna$seurat_clusters, 
          from=seq(0, 31), 
          to=new.cluster.ids)

# Prepare data for integration
seurat_rna <- NormalizeData(seurat_rna)
seurat_rna <- ScaleData(seurat_rna, features = rownames(seurat_rna))
seurat_rna <- FindVariableFeatures(seurat_rna)
seurat_rna <- RunPCA(seurat_rna)
seurat_rna <- RunUMAP(seurat_rna, dims = 1:30)

# Label trasnfer with ArchR
proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seurat_rna,
    addToArrow = TRUE,
    groupRNA = "annotated_clusters",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
)

# Get cell identity by cluster and name it by most frequent cell type in cluster
cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
proj$Clusters2 <- mapvalues(proj$Clusters, 
          from=rownames(cM), 
          to=labelNew)

# Assign OPCs clusters
proj$Clusters2[which(proj$predictedGroup_Un == "OPCs")] <- "OPCs"
pal <- paletteDiscrete(values = seurat_rna$seurat_clusters)

# Plot annotated UMAP
p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony")

# Differential Gene activity
# OPCs vs pre-OPCs
markers <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = "OPCs",
    bgdGroups = "Glioblast/Pre-OPC"
)

# One vs all
markers_all <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markers, cutOff = "FDR <= 0.1")
markerList_all <- getMarkers(markers_all, cutOff = "FDR <= 0.1")

preopcs_markers <- head(markerList[[1]]$name, 30)
opcs_markers <- markerList[[1]]

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    useGroups = c("OPCs", "Glioblast/Pre-OPC"),
    geneSymbol = opcs_markers, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30,
    sizes = c(2, 2, 2, 2)

)

# Find peaks with MACS2
proj<- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters2")
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(proj)
proj <- addPeakMatrix(proj)

# Find diffentiable accessible peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = "OPCs",
    bgdGroups = "Glioblast/Pre-OPC"
)
markersPeaks_all <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList_peaks <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.0")
markerList_peaks_down <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC <= 0.0")

# ChromVAR TF motif deviations
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# Get differential deviations OPCs vs pre-OPCs
markersMotifs <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "MotifMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = "OPCs",
    bgdGroups = "Glioblast/Pre-OPC",
    useSeqnames = "z"
)
proj_opcs_preopcs <- proj[which(proj$Clusters2 == "OPCs" | proj$Clusters2 == "Glioblast/Pre-OPC"), ]

# Save Archr Project
saveRDS(file='/data/proj/GCB_FBP/human_dev/archr_proj_best3.rds', proj)
proj <- readRDS(file='/data/proj/GCB_FBP/human_dev/archr_proj_best3.rds')

# Plots
# Customize UMAP plot
umap_coord <- proj@embeddings[[2]][[1]]
colnames(umap_coord) <- c("UMAP_1", "UMAP_2")
umap_coord$samples <- proj$Sample
umap_coord$clusters <- proj$Clusters2
umap_coord_2 <- cbind(umap_coord, celltype.predictions) 
proj_path <- "/data/proj/GCB_FBP/human_dev/exploratory/hca_dev_best3/Plots/"

## Match colors with scRNA 
color_df <- data.frame(seurat_clusters = seurat_rna$seurat_clusters, clusters_names = seurat_rna$annotated_clusters)
color_df <- arrange(color_df[!duplicated(color_df$seurat_clusters), ], seurat_clusters)
color_df$color <- sidecols
color_df <- color_df[color_df$clusters_names %in% as.factor(unique(proj$Clusters2)), ]
color_df$clusters_names <- as.character(color_df$clusters_names)
color_df$clusters_names <- factor(color_df$clusters_names)
levels(color_df$clusters_names) <- levels(factor(proj$Clusters2))
color_df <- arrange(color_df, clusters_names)
pal <- c("#7C7328","#C83ADB","#8CA35F","#762B29","#D840A9","#434889","#8B5326","#47A78A","#49B532","#E14121","#3A712C","#DC7756","#A33424","#D176AC","#7D317A","#DF7928","#6932DE","#BE65D5","#D84076","#499EBD","#C59131","#532796","#C09064","#8E3155","#A5A330","#A583CD","#628DD0","#73A330","#6365DA","#4AAD61","#CE787F","#DF3E4D")

sidecols <- iwanthue(length(unique(networkAnnotableFile$seurat_clusters)), 0, 360, 36, 180, 13, 73)
sidecols <- sidecols[as.factor(networkAnnotableFile$seurat_clusters)]
names(new.cluster.ids) <- levels(umap_coord$clusters)

umap_coord$ab_clusters <- mapvalues(umap_coord$clusters, 
          from=unique(umap_coord$clusters), 
          to=new.cluster.ids)

pal <- paletteDiscrete(values = umap_coord$ab_clusters)
pal <- color_df$color
p <- ggplot(umap_coord, aes(UMAP_1, UMAP_2, color = factor(clusters))) +
    geom_point(size=0.1, alpha=0.9) + #, show.legend = FALSE) +
    theme_cowplot(12) +
    scale_color_manual(values=pal) 

pdf("Custom-UMAP-clusters.pdf", 10, 10)
ggplot(umap_coord, aes(UMAP_1, UMAP_2, color = factor(ab_clusters))) +
    geom_point(size=1, alpha=0.7) +
    theme_cowplot(12) +
    scale_color_manual(values=pal) 
dev.off()

pal <- paletteDiscrete(values = umap_coord$samples)
p <- ggplot(umap_coord, aes(UMAP_1, UMAP_2, color = factor(samples))) +
    geom_point(size=0.1, alpha=0.9) + #, show.legend = FALSE) +
    theme_cowplot(12) +
    scale_color_manual(values=pal) 
plotPDF(p, name = "Plot-UMAP-Sample-Age.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

p <- ggplot(umap_coord_2, aes(UMAP_1, UMAP_2, color = prediction.score.OPCs)) +
    geom_point(size=0.1, alpha=0.9) + #, show.legend = FALSE) +
    theme_cowplot(12) 
pdf(file = paste0(proj_path, "Plot-UMAP-OPCs-Score.pdf"))
p
dev.off()

p <- ggplot(umap_coord_2, aes(UMAP_1, UMAP_2, color = prediction.score.Glioblast.Pre.OPC)) +
    geom_point(size=0.1, alpha=0.9) + #, show.legend = FALSE) +
    theme_cowplot(12) 
pdf(file = paste0(proj_path, "Plot-UMAP-PreOPCs-Score.pdf"))
p
dev.off()

# Proportions barplot
prop <- data.frame(clusters = proj$Clusters2, sampleid = proj$Sample)
p <- ggplot(prop, aes(clusters, fill = sampleid)) +
    geom_bar(position="fill") +
    theme_cowplot(12) +
    scale_fill_manual(values=pal) +
    coord_flip()
plotPDF(p, name = "Plot-Cluster-proportion.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

# Violin plots RNA
seurat_rna_opcs <- subset(seurat_rna, subset = annotated_clusters == "OPCs" | annotated_clusters == "Glioblast/Pre-OPC")
pdf(file = paste0(proj_path, "primed_genes_expression.pdf"), 10, 5)
features = c("LUZP2", "PDGFRA", "ETV1")
VlnPlot(seurat_rna_opcs, features = features, group.by = "annotated_clusters")
dev.off()

# Plot tracks
# Compare gene activity and track plots from genes from scRNA analysis
gene_uni <- getFeatures(proj)
opcs <- c("CXXC4","NR2F2","MEST","OTX2","GAD2","MIR99AHG","WLS","AL139246.5","OLIG2","PAX3","PBX3","IRX1",
          "DMBX1","CRNDE","OTX2-AS1","RMST","TOX3","IRX2","GAD1","LMO1","ANKRD6","MGST1","DLX1","LINC01896",
          "ZIC3","PLP1","LINC00237","MAB21L1","SOX14","OLIG1","SOX1-OT","LHX5-AS1","ZIC1","DLX2","EDNRB",
          "CCND1","FZD10-AS1","LRRN3","IRX5","ID3","GSX2","ASCL1","IRX3","FZD10","ZBTB20","SALL3","SEMA6A",
          "AC233296.1","PDZRN3","SOX2-OT","LRP2","LRRC4C","HELT","SPATS2L","SLITRK2","KIZ","BCHE","CEP112",
          "LITAF","LINC01210","SP9","NKX2-1","C5orf38","ZIC4","SLIT2","CHL1","PROX1","LHX5","AC067773.1",
          "HS3ST3A1","SP5","LGR4","NTRK2","EFNA5","EPB41L4B","GSX1","PTCH1","CT75","MECOM","SHD","ASIC4",
          "GATA3","TPPP3","AC087477.2","SIX3","SCRG1","EGFEM1P","HOXB2","PGM5","CDKN1C","NES","AC090152.1",
          "ST8SIA5","STK33","AK4","DIPK1C","EGFR","SH3BP4","LHX1", "GNG5")
opcs <- opcs[which(opcs %in% gene_uni)]
opcs_degs <- opcs[which(opcs %in% opcs_markers$name)]
opcs_degs <- opcs_markers[which(opcs_markers$name %in% opcs), ]

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    useGroups = c("OPCs", "Glioblast/Pre-OPC"),
    geneSymbol = opcs, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30
)
plotPDF(p, name = "Plot-Tracks-OPCs-rna", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    geneSymbol = opcs, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30
)
plotPDF(p, name = "Plot-Tracks-AllClusters-OPCs-rna", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = opcs, 
    embedding = "UMAPHarmony",
    imputeWeights = NULL
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-OPCs-rnas-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = opcs, 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(proj)

)
plotPDF(plotList = p, 
    name = "Plot-UMAP-OPCs-rnas-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

opcs_restricted <- c("S100B","APOD","SCRG1","PDGFRA","DBI","OLIG1","OLIG2","COL9A3","BAMBI","C2orf27A", 
                    "SERPINE2","LIMA1","RGCC","ITM2C","COL9A1","BCAN","PLPPR1","MT3","PMP2","SOX10",
                    "KLRC2","NKX2-2","PTN","RAMP1","EDNRB","ASB3","LRRC4C","PLLP","LUZP2","NCALD",
                    "LINC01896","FABP7","CMTM5","ASIC4","PLPP4","RIT2","IDH1","PLAT","SIRT2","SCD5",
                    "ARL4A","CA10","GPM6A","NXPH1","BCHE","NKAIN4","TMEM121","C1orf61","KCNIP1","RAB31",
                    "PID1","LINC00643","SHISA4","RHOC","RPRM","PTPRZ1","GPR17","CNP","OMG","TSPAN7",
                    "DNER","AC091138.1","GRB14","PDE4B","SOX6","CNTN1","SNTG1","CSPG5","PPP2R2B","B2M",
                    "COL20A1","TTYH1","COTL1","FAM89A","LSAMP","DOCK10","NENF","BRINP1","ETV1","ARL2BP",
                    "CADM2","APLP2","LAMA4","TRAF4","CNMD","LINC00511","CTHRC1","C1QL1","TAOK3","TMEM9B",
                    "MT2A","TIMP4","NTRK3","SLC22A17","BAALC","MAGEH1","FIP1L1","FBXO7","KANK1","PCDH9")
res_degs <- opcs_markers[which(opcs_markers$name %in% opcs_restricted), ]

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    useGroups = c("OPCs", "Glioblast/Pre-OPC"),
    geneSymbol = res_degs$name, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30,
    sizes = c(2, 2, 2, 2),
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC > 0", returnGR = TRUE)["OPCs"] 
)
plotPDF(p, name = "Plot-Tracks-DEGs-OPCs-restricted-rna", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

pal <- paletteContinuous(set = "greyMagma", n=256, reverse=FALSE)
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = res_degs$name, 
    embedding = "UMAPHarmony",
    imputeWeights = NULL,
    pal = pal
)
plotPDF(p, name = "Plot-UMAP-DEGs-OPCs-restricted-rnas-Marker-Genes-WO-Imputation.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

tf_opcs_ol <- c("OLIG1","OLIG2","SOX10","NKX2-2","LUZP2","NCALD","KCNIP1","SOX6","ETV1","TRAF4","FIP1L1","PDLIM5",
                "TPI1","ZCCHC17","JUNB","UGP2","ETV5","NR0B1","ETV4","ENO1","SUCLG1","EPAS1","FOSB","TSC22D4",
                "XBP1","HMX1","EGR1","CERS4","ETS1","TCEAL2","SALL3","MITF","DAB2","FOS","NME1","HEY1",
                "CREB5","CANX","PRNP","RAB2A","NPDC1","CKMT1B","HOXB2","RAB7A","PPP5C","CD59","ZNF77","KLF4",
                "CREM","CFL2","RBM42","HMGA1","AHR","KLF6","UBB","KIF22","HOXA3","ZNF580","ERF","PTPMT1",
                "HMGN3","JDP2","DNMT1","BAX","HIVEP1","MYNN","SOX8","CREB3L1","IKZF2","NFATC1","MTHFD1","TCF7L2",
                "ZBTB47","ZNF607","ZNF44","CLK1","NR6A1","ACO1","NR3C2","ZNF625","MAFF","ZNF181","TFDP1","BATF3",
                "NUCB1","SFT2D1","NUP133","NMI","FOXO3","NFE2L2","MAP4K2","KLF5","RFX4","EGR2","TRIM24","ZBTB39",
                "ZNF649","GRHPR","PPARGC1A","GAR1")
tf_degs <- opcs_markers[which(opcs_markers$name %in% tf_opcs_ol), ]

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    useGroups = c("OPCs", "Glioblast/Pre-OPC"),
    geneSymbol = tf_opcs_ol, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30,
    sizes = c(2, 2, 2, 2)

)
plotPDF(p, name = "Plot-Tracks-OPCs-TF-rna", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = tf_opcs_ol, 
    embedding = "UMAPHarmony",
    imputeWeights = NULL,
    pal = pal
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-OPCs-TF-rnas-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

tf_opcs_all <- c("NR2F2","OTX2","OLIG2","PAX3","PBX3","IRX1","DMBX1","IRX2","DLX1","ZIC3","SOX14","OLIG1","ZIC1",
                 "DLX2","IRX5","GSX2","ASCL1","IRX3","ZBTB20","SALL3","HELT","SP9","NKX2-1","ZIC4","PROX1","LHX5",
                 "SP5","GSX1","MECOM","GATA3","SIX3","HOXB2","LHX1","DLX5","PDLIM5","NKX2-2","HOXB3","ID1","PAX7",
                 "OLIG3","TCF7L2","HOXA3","SOX2","MSI2","RFX4","FOXP2","LMO2","FEZF1","TFDP2","GBX2","DLX6","ZIC2",
                 "UNCX","HOXB4","SRY","HOXA2","ZNF503","VAX1","ESRRG","OTP","LHX8","MEIS1","HOXC4","ATOH1","HOXD3",
                 "EN2","ZNF385A","FOXO1","ZFHX3","HOXA4","ZIC5","SOX13","ARID5B","PRDM12","ZEB1","CTNNB1","SOX8","ZNF521",
                 "IRX4","NKX2-3","ETFB","EN1","HES5","HMGB2","TAL1","LUZP2","HOXA5","FANK1","HOXB5","FOXB1","SOX21",
                 "PITX2","SOX10","SOX6","NKX2-4","FOXA1","EGR1","POU3F2","NKX6-1","HOXD4")

tf_all_degs <- opcs_markers[which(opcs_markers$name %in% tf_opcs_all), ]

p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters2",
    useGroups = c("OPCs", "Glioblast/Pre-OPC"),
    geneSymbol = tf_opcs_all, 
    minCells = 0,
    upstream = 50000,
    downstream = 50000,
    plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
    scCellsMax = 30
)

tf_opcs_all <- tf_opcs_all[which(tf_opcs_all %in% gene_uni)]
plotPDF(p, name = "Plot-Tracks-OPCs-TFall-rna", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = tf_opcs_all, 
    embedding = "UMAPHarmony",
    imputeWeights = NULL
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-OPCs-TFall-rnas-Marker-Genes-WO-Imputation.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 5, height = 5)

# ChromVARdeviations plot
# Smoothing
proj_opcs_preopcs <- addImputeWeights(proj_opcs_preopcs)

p <- plotGroups(ArchRProj = proj_opcs_preopcs, 
  groupBy = "Clusters2", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj_opcs_preopcs)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})

pdf(paste0(proj_path,"chromvar_deviations.pdf"), 100, 3)
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
dev.off()


