############################################################################################################
# Reads in important libraries for downstream analysis
library(tidyverse)

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

library(SoupX)

library(pcaMethods)
library(hdf5r)
library(harmony)
library(rliger)
library(leiden)

library(monocle3)

library(URD)
library(Matrix)
library(rgl)

library(velociraptor)

library(ShinyCell)

library(ggplot2)
library(ggpointdensity)
library(viridis)
library(future)
library(sctransform)
library(Rcpp)
library(cowplot)
library(patchwork)
library(dittoSeq)
library(stringr)
library(knitr)
library(clustree)
library(stats)
library(scater)
library(gtools)
library(magrittr)
library(grid)
library(reticulate)
library(dplyr)
#https://distill.pub/2016/misread-tsne/
#https://pair-code.github.io/understanding-umap/
############################################################################################################





############################################################################################################
FindOutlierThershold <- function(Meta.Data)
{
  low.outlier.thershold <- median(Meta.Data) - 3*mad(Meta.Data)
  high.outlier.thershold <- median(Meta.Data) + 3*mad(Meta.Data)
  
  return(list(low.outlier.thershold, high.outlier.thershold))
}
############################################################################################################





############################################################################################################
GenomeSpecificCellCycleScoring <- function(Seurat.obj, Genome)
{
  # Pull standard cell cycle genes for cell cycle scoring
  human.s.genes <- cc.genes.updated.2019$s.genes
  human.g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  # Had to append "hg19-" to the beginning of all genes due to dual index reference genome used
  dual.hg19.s.genes <- paste("hg19-", human.s.genes, sep="")
  dual.hg19.g2m.genes <- paste("hg19-", human.g2m.genes, sep="")
  
  # Had to append "GRCh38-" to the beginning of all genes due to dual index reference genome used
  dual.GRCh38.s.genes <- paste("GRCh38-", human.s.genes, sep="")
  dual.GRCh38.g2m.genes <- paste("GRCh38-", human.g2m.genes, sep="")
  
  # Changed genes to title case for mice genes (Human genes are annotated in all CAPS)
  mm10.s.genes <- str_to_title(cc.genes.updated.2019$s.genes)
  mm10.g2m.genes <- str_to_title(cc.genes.updated.2019$g2m.genes)
  
  # Appended "GRCh38-" to the beginning of all genes due to dual index reference genome used
  dual.mm10.s.genes <- paste("mm10---", mm10.s.genes, sep="")
  dual.mm10.g2m.genes <- paste("mm10---", mm10.g2m.genes, sep="")
  
  
  if (Genome == "GRCh38" | Genome == "hg19"){
    Seurat.obj <- CellCycleScoring(Seurat.obj, s.features = human.s.genes, g2m.features = human.g2m.genes)
  }
  else if (Genome == "mm10"){
    Seurat.obj <- CellCycleScoring(Seurat.obj, s.features = mm10.s.genes, g2m.features = mm10.g2m.genes)
  }
  else if (Genome == "Dualhg19"){
    Seurat.obj <- CellCycleScoring(Seurat.obj, s.features = dual.hg19.s.genes, g2m.features = dual.hg19.g2m.genes)
  }
  else if (Genome == "DualGRCh38"){
    Seurat.obj <- CellCycleScoring(Seurat.obj, s.features = dual.GRCh38.s.genes, g2m.features = dual.GRCh38.g2m.genes)
  }
  else if (Genome == "Dualmm10"){
    Seurat.obj <- CellCycleScoring(Seurat.obj, s.features = dual.mm10.s.genes, g2m.features = dual.mm10.g2m.genes)
  }
}
############################################################################################################





############################################################################################################
Filter.GEMs <- function(Seurat.obj, Project.Path, GEM.Path, Figure.Total, Analysis.Name, Genome)
{
  # Reading in sample specific GEM classification data outputted from CellRanger count when aligned to a Dual Reference Genome
  GEM.classification.data <- read.csv(paste(GEM.Path, "/analysis/gem_classification.csv", sep = ""), sep = ',', row.names = 1)
  Seurat.obj <- AddMetaData(Seurat.obj, GEM.classification.data, col.name = paste("GEM.classification." , colnames(GEM.classification.data), sep = ""))
  
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_GEM_Filtering/", sep = ""))
  
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_GEM_Filtering/", Analysis.Name, "_GEM_classification_graphs.pdf", sep = ""), width = 12, height = 12)
  
  if("hg19" %in% Genome)
  {
    print(ggplot(GEM.classification.data, aes(x=hg19, y=mm10, col=call)) + geom_point(alpha=0.5, shape=1) + scale_color_manual(values = c("black", "#E69F00", "#56B4E9")))
    print(ggplot(subset(GEM.classification.data, call == "hg19"), aes(x=hg19, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="black"))
    print(ggplot(subset(GEM.classification.data, call == "mm10"), aes(x=hg19, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="#E69F00"))
    print(ggplot(subset(GEM.classification.data, call == "Multiplet"), aes(x=hg19, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="#56B4E9"))
    print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "GEM.classification.call") + theme_bw() + theme(aspect.ratio = 1))
    
    human.cells <- subset(Seurat.obj, subset = GEM.classification.call == Genome)
    print(FeatureScatter(object = human.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "red") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Human cells"))
    
    mouse.cells <- subset(Seurat.obj, subset = GEM.classification.call == "mm10")
    print(FeatureScatter(object = mouse.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "green") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Mouse cells"))
    
    multiplet.cells <- subset(Seurat.obj, subset = GEM.classification.call == "Multiplet")
    print(FeatureScatter(object = multiplet.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "blue") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Multiplet cells"))
    
    dev.off()
  }
  else if("GRCh38" %in% Genome)
  {
    print(ggplot(GEM.classification.data, aes(x=GRCh38, y=mm10, col=call)) + geom_point(alpha=0.5, shape=1) + scale_color_manual(values = c("black", "#E69F00", "#56B4E9")))
    print(ggplot(subset(GEM.classification.data, call == "GRCh38"), aes(x=GRCh38, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="black"))
    print(ggplot(subset(GEM.classification.data, call == "mm10"), aes(x=GRCh38, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="#E69F00"))
    print(ggplot(subset(GEM.classification.data, call == "Multiplet"), aes(x=GRCh38, y=mm10, col=call)) + geom_point(alpha=1, shape=1, color="#56B4E9"))
    print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "GEM.classification.call") + theme_bw() + theme(aspect.ratio = 1))
    
    human.cells <- subset(Seurat.obj, subset = GEM.classification.call == Genome)
    print(FeatureScatter(object = human.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "red") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Human cells"))
    
    mouse.cells <- subset(Seurat.obj, subset = GEM.classification.call == "mm10")
    print(FeatureScatter(object = mouse.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "green") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Mouse cells"))
    
    multiplet.cells <- subset(Seurat.obj, subset = GEM.classification.call == "Multiplet")
    print(FeatureScatter(object = multiplet.cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = "blue") + theme_bw() + theme(aspect.ratio = 1) + ggtitle("Multiplet cells"))
    
    dev.off()
  }
  
  
  write.csv(summary(GEM.classification.data), file=paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_GEM_Filtering/", Analysis.Name, "_GEM_classification_stat_summary.csv", sep = ""))
  
  # Removing cells (GEMs) that CellRanger identified as mouse (majority of sequences aligned to mm10)
  Seurat.obj <- subset(Seurat.obj, subset = GEM.classification.call == Genome)
  
  return(Seurat.obj)
}
############################################################################################################





############################################################################################################
Generate.QC.Plots <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, File.Name)
{
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", File.Name, ".pdf", sep = ""), width = length(unique(Seurat.obj$ID))*12, height = 12)
  
  print(ggplot() + geom_density(data=Seurat.obj@meta.data, mapping=aes(color=ID, x=nCount_RNA, fill=ID), alpha=0.2) + scale_x_log10() + theme_classic() + ylab("Cell Density") + theme(aspect.ratio = 1))
  
  print(ggplot() + geom_density(data=Seurat.obj@meta.data, mapping=aes(color=ID, x=nFeature_RNA, fill=ID), alpha=0.2) + scale_x_log10() + theme_classic() + theme(aspect.ratio = 1))
  
  print(ggplot() + geom_point(data=Seurat.obj@meta.data, mapping=aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mito)) + scale_color_gradient(low="grey90", high="black") + theme_bw() + theme(aspect.ratio = 1))
  
  print(ggplot() + geom_density(data=Seurat.obj@meta.data, mapping=aes(color=ID, x=percent.mito, fill=ID), alpha=0.2) + scale_x_log10() + theme_classic() + theme(aspect.ratio = 1))
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  print(ggplot() + geom_density(data=Seurat.obj@meta.data, mapping=aes(x=log10GenesPerUMI, color=ID, fill=ID), alpha=0.2) + theme_classic() + theme(aspect.ratio = 1))
  
  print(VlnPlot(object = Seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), group.by = "orig.ident", ncol = 3))
  print(VlnPlot(object = Seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), group.by = "ID", ncol = 3))
  
  print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "ID") + theme_bw() + theme(aspect.ratio = 1))
  print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mito", group.by = "ID") + geom_bin_2d(bins = 300) + scale_fill_continuous(type = "viridis") + theme_bw() + theme(aspect.ratio = 1))

    print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "ID") + theme_bw() + theme(aspect.ratio = 1))
    print(FeatureScatter(object = Seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "ID") + geom_bin_2d(bins = 300) + scale_fill_continuous(type = "viridis") + theme_bw() + theme(aspect.ratio = 1))

    print(ggplot() + geom_bar(data=Seurat.obj@meta.data, mapping=aes(x=ID, fill=ID)) + theme_classic() + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) + theme(plot.title=element_text(hjust=0.5, face="bold")) + ggtitle("Number of Cells") + geom_label(aes(x=Seurat.obj$ID, y=length(Seurat.obj$ID), label=length(Seurat.obj$ID))))
  
  dev.off()
}
############################################################################################################





############################################################################################################
Single.Sample <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, Genome, Regress.Cell.Cycle, PCA.Feature.List, Resolution.List)
{
  Seurat.obj <- NormalizeData(Seurat.obj)
  Seurat.obj <- FindVariableFeatures(Seurat.obj, selection.method = "vst", nfeatures = 3000)
  
  Seurat.obj <- GenomeSpecificCellCycleScoring(Seurat.obj, Genome)
  
  if (Regress.Cell.Cycle == "NO"){
    Seurat.obj <- ScaleData(Seurat.obj, features = rownames(Seurat.obj))
  }
  else if (Regress.Cell.Cycle == "YES"){
    Seurat.obj <- ScaleData(Seurat.obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.obj))
  }
  else if (Regress.Cell.Cycle == "DIFF"){
    Seurat.obj$CC.Difference <- Seurat.obj$S.Score - Seurat.obj$G2M.Score
    Seurat.obj <- ScaleData(Seurat.obj, vars.to.regress = "CC.Difference", features = rownames(Seurat.obj))
  }
  
  
  
  if (missing(PCA.Feature.List)){
    Seurat.obj <- RunPCA(Seurat.obj, features = VariableFeatures(Seurat.obj), verbose = FALSE)
  }
  else{
    Seurat.obj <- RunPCA(Seurat.obj, features = PCA.Feature.List, verbose = FALSE)
  }
  
  Seurat.obj <- RunUMAP(Seurat.obj, dims = 1:30)
  Seurat.obj <- FindNeighbors(Seurat.obj, dims = 1:30)
  
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Single.Sample/", sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Single.Sample/", Analysis.Name, "_all_cluster_resolutions.pdf", sep = ""), width = 12, height = 12)
  
  
  if (missing(Resolution.List)){
    for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.obj <- FindClusters(Seurat.obj, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.obj, label = TRUE) + ggtitle(paste("Resolution: ", res)))
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Single.Sample/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.obj, prefix = "RNA_snn_res."))
    print(clustree(Seurat.obj, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) > 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.obj <- FindClusters(Seurat.obj, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.obj, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Single.Sample/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.obj, prefix = "RNA_snn_res."))
    print(clustree(Seurat.obj, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) == 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.obj <- FindClusters(Seurat.obj, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.obj, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
  }
  
  return(Seurat.obj)
}
############################################################################################################





############################################################################################################
Simple.Merge <- function(Seurat.List, Project.Path, Figure.Total, Analysis.Name, Genome, Regress.Cell.Cycle, PCA.Feature.List, Resolution.List)
{
  
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x) {
    x <- NormalizeData(x)
  })
  Seurat.combined <- merge(Seurat.List[[1]], Seurat.List[c(-1)])
  Seurat.combined <- FindVariableFeatures(Seurat.combined, selection.method = "vst", nfeatures = 3000)
  
  Seurat.combined <- GenomeSpecificCellCycleScoring(Seurat.combined, Genome)
  
  if (Regress.Cell.Cycle == "NO"){
    Seurat.combined <- ScaleData(Seurat.combined, features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "YES"){
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "DIFF"){
    Seurat.combined$CC.Difference <- Seurat.combined$S.Score - Seurat.combined$G2M.Score
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = "CC.Difference", features = rownames(Seurat.combined))
  }
  
  
  
  if (missing(PCA.Feature.List)){
    Seurat.combined <- RunPCA(Seurat.combined, features = VariableFeatures(Seurat.combined), verbose = FALSE)
  }
  else{
    Seurat.combined <- RunPCA(Seurat.combined, features = PCA.Feature.List, verbose = FALSE)
  }
  
  Seurat.combined <- RunUMAP(Seurat.combined, dims = 1:30)
  Seurat.combined <- FindNeighbors(Seurat.combined, dims = 1:30)
  
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Simple.Merge/", sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Simple.Merge/", Analysis.Name, "_all_cluster_resolutions.pdf", sep = ""), width = 12, height = 12)
  
  if (missing(Resolution.List)){
    for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)))
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Simple.Merge/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) > 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Simple.Merge/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) == 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
  }
  
  return(Seurat.combined)
}
############################################################################################################





############################################################################################################
Seurat.Integration <- function(Seurat.List, Project.Path, Figure.Total, Analysis.Name, Genome, Regress.Cell.Cycle, PCA.Feature.List, Resolution.List, Big.Data, Reference.List)
{
  cat("Beginning Seurat Integration pipeline for: ", Analysis.Name, "\n")
  
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x)
  {
    cat("Conducting SCTransform for sample:", unique(x$ID), "\n")
    x <- SCTransform(x, variable.features.n = 3000, verbose = TRUE, return.only.var.genes = FALSE)
  })
  
  cat("Selecting Intgration Features\n")
  features <- SelectIntegrationFeatures(object.list = Seurat.List, nfeatures = 3000)
  Seurat.List <- PrepSCTIntegration(object.list = Seurat.List, anchor.features = features)
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x)
  {
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
    
  if (Big.Data == TRUE)
  {
    if (missing(Reference.List))
    {
      anchors <- FindIntegrationAnchors(object.list = Seurat.List, normalization.method = "SCT", anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
      Seurat.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)
    }
    else
    {
      anchors <- FindIntegrationAnchors(object.list = Seurat.List, reference = Reference.List, normalization.method = "SCT", anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 20)
      Seurat.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:50)
    }
  }
  else
  {
    if (missing(Reference.List))
    {
      anchors <- FindIntegrationAnchors(object.list = Seurat.List, normalization.method = "SCT", anchor.features = features)
      Seurat.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    }
    else
    {
      anchors <- FindIntegrationAnchors(object.list = Seurat.List, reference = Reference.List, normalization.method = "SCT", anchor.features = features)
      Seurat.combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    }
  }
  
  Seurat.combined <- GenomeSpecificCellCycleScoring(Seurat.combined, Genome)
  
  if (Regress.Cell.Cycle == "NO"){
    Seurat.combined <- ScaleData(Seurat.combined, features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "YES"){
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "DIFF"){
    Seurat.combined$CC.Difference <- Seurat.combined$S.Score - Seurat.combined$G2M.Score
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = "CC.Difference", features = rownames(Seurat.combined))
  }
  
  
  
  if (missing(PCA.Feature.List)){
    Seurat.combined <- RunPCA(Seurat.combined, features = VariableFeatures(Seurat.combined), verbose = FALSE)
  }
  else{
    Seurat.combined <- RunPCA(Seurat.combined, features = PCA.Feature.List, verbose = FALSE)
  }
  
  Seurat.combined <- RunUMAP(Seurat.combined, reduction = "pca", dims = 1:30)
  Seurat.combined <- FindNeighbors(Seurat.combined, dims = 1:30)
  
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Seurat.Integration/", sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Seurat.Integration/", Analysis.Name, "_all_cluster_resolutions.pdf", sep = ""), width = 12, height = 12)
  
  if (missing(Resolution.List)){
    for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)))
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Seurat.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) > 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Seurat.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) == 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
  }
  
  return(Seurat.combined)
}
############################################################################################################





############################################################################################################
Harmony.Integration <- function(Seurat.List, Project.Path, Figure.Total, Analysis.Name, Genome, Regress.Cell.Cycle, PCA.Feature.List, Resolution.List)
{
  cat("Beginning Harmony Integration pipeline for: ", Analysis.Name, "\n")
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x) {
    cat("Normalizing Data for sample: ", unique(x$ID), "\n")
    x <- NormalizeData(x)
  })
  cat("Merging all Seurat Objects into a single Seurat Object\n")
  Seurat.combined <- merge(Seurat.List[[1]], Seurat.List[c(-1)])
  cat("Identifying highly variable features\n")
  Seurat.combined <- FindVariableFeatures(Seurat.combined, selection.method = "vst", nfeatures = 3000)
  
  cat("Assigning Cell-Cycle Scores\n")
  Seurat.combined <- GenomeSpecificCellCycleScoring(Seurat.combined, Genome)
  
  if (Regress.Cell.Cycle == "NO"){
    cat("Scaling the data\n")
    Seurat.combined <- ScaleData(Seurat.combined, features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "YES"){
    cat("Regressing out Cell Cycle\n")
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.combined))
  }
  else if (Regress.Cell.Cycle == "DIFF"){
    cat("Regressing out the difference between the G2M and S phase scores\n")
    Seurat.combined$CC.Difference <- Seurat.combined$S.Score - Seurat.combined$G2M.Score
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = "CC.Difference", features = rownames(Seurat.combined))
  }
  
  
  
  if (missing(PCA.Feature.List)){
    cat("Performing linear dimensional reduction on highly variable genes\n")
    Seurat.combined <- RunPCA(Seurat.combined, features = VariableFeatures(Seurat.combined), verbose = FALSE)
  }
  else{
    cat("Performing linear dimensional reduction on custom gene list\n")
    Seurat.combined <- RunPCA(Seurat.combined, features = PCA.Feature.List, verbose = FALSE)
  }
  
  cat("Performing Harmony Integration\n")
  Seurat.combined <- RunHarmony(Seurat.combined, group.by.vars = "ID", plot_convergence = TRUE)
  
  cat("Performing Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique\n")
  Seurat.combined <- RunUMAP(Seurat.combined, reduction = "harmony", dims = 1:30)
  cat("Computing the k.param nearest neighbors for a given dataset\n")
  Seurat.combined <- FindNeighbors(Seurat.combined, reduction = "harmony", dims = 1:30)
  
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Harmony.Integration/", sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Harmony.Integration/", Analysis.Name, "_all_cluster_resolutions.pdf", sep = ""), width = 12, height = 12)
  
  
  if (missing(Resolution.List)){
    for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)))
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Harmony.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) > 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Harmony.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) == 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
  }
  
  cat("Finishing Harmony Integration pipeline for: ", Analysis.Name, "\n")
  return(Seurat.combined)
}
############################################################################################################





############################################################################################################
Liger.Integration <- function(Seurat.List, Project.Path, Figure.Total, Analysis.Name, Genome, Regress.Cell.Cycle, Resolution.List)
{
  
  Seurat.List <- lapply(X = Seurat.List, FUN = function(x) {
    x <- NormalizeData(x)
  })
  Seurat.combined <- merge(Seurat.List[[1]], Seurat.List[c(-1)])
  Seurat.combined <- FindVariableFeatures(Seurat.combined, selection.method = "vst", nfeatures = 3000)
  
  Seurat.combined <- GenomeSpecificCellCycleScoring(Seurat.combined, Genome)
  
  if (Regress.Cell.Cycle == "NO"){
    Seurat.combined <- ScaleData(Seurat.combined, features = rownames(Seurat.combined), split.by = "ID", do.center = FALSE)
  }
  else if (Regress.Cell.Cycle == "YES"){
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.combined), split.by = "ID", do.center = FALSE)
  }
  else if (Regress.Cell.Cycle == "DIFF"){
    Seurat.combined$CC.Difference <- Seurat.combined$S.Score - Seurat.combined$G2M.Score
    Seurat.combined <- ScaleData(Seurat.combined, vars.to.regress = "CC.Difference", features = rownames(Seurat.combined), split.by = "ID", do.center = FALSE)
  }
  
  
  
  Seurat.combined <- RunOptimizeALS(Seurat.combined, k = 20, lambda = 5, split.by = "ID")
  Seurat.combined <- RunQuantileNorm(Seurat.combined, split.by = "ID")
  # You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
  # `RunQuantileNorm` according to your needs
  Seurat.combined <- FindNeighbors(Seurat.combined, reduction = "iNMF", dims = 1:20)
  # Dimensional reduction and plotting
  Seurat.combined <- RunUMAP(Seurat.combined, dims = 1:ncol(Seurat.combined[["iNMF"]]), reduction = "iNMF")
  
  dir.create(paste(Project.Path, "/", Analysis.Name, sep = ""))
  dir.create(paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Liger.Integration/", sep = ""))
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Liger.Integration/", Analysis.Name, "_cluster_resolution_", res, ".pdf", sep = ""), width = 12, height = 12)
  
  
  if (missing(Resolution.List)){
    for (res in c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)))
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Liger.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) > 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
    
    pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_Liger.Integration/", Analysis.Name, "_cluster_tree.pdf", sep = ""), width = 60, height = 60)
    print(clustree(Seurat.combined, prefix = "RNA_snn_res."))
    print(clustree(Seurat.combined, prefix = "RNA_snn_res.", node_colour = "sc3_stability"))
    dev.off()
  }
  else if (length(Resolution.List) == 1){
    for (res in Resolution.List)
    {
      message("FindClusters is being calculated for resolution:", res)
      Seurat.combined <- FindClusters(Seurat.combined, resolution = res, algorithm = 4, method = "igraph")
      print(DimPlot(Seurat.combined, label = TRUE) + ggtitle(paste("Resolution: ", res)) )
    }
    dev.off()
  }
  
  return(Seurat.combined)
}
############################################################################################################





############################################################################################################
Generate.Differential.Data <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, Integration.Method)
{
  all.markers <- FindAllMarkers(Seurat.obj, assay = "RNA", logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1, return.thresh = 0.01)
  write.table(all.markers, file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_All_markers.tsv", sep = ""), sep = "\t", quote = FALSE)
  
  all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
  write.table(top10, file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_All_markers_Top10.tsv", sep = ""), sep = "\t", quote = FALSE)
  
  all.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC) -> top2
  write.table(top2, file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_All_markers_Top2.tsv", sep = ""), sep = "\t", quote = FALSE)
  
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_All_markers_Top10_Heatmap.pdf", sep = ""), width = 12, height = 12)
  print(DoHeatmap(Seurat.obj, features = top10$gene))
  dev.off()
}


############################################################################################################





############################################################################################################
Generate.Feature.Plots <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, Integration.Method, Cell.Type.DF, Genome)
{
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_feature_plots.pdf", sep = ""), width = 6, height = 6)
  
  for (i in 1:ncol(Cell.Type.DF))
  {
    cell.type <- colnames(Cell.Type.DF)[i]
    gene.markers <- as.character(Cell.Type.DF[, i][Cell.Type.DF[, i] != ""])
    
    if (Genome == "GRCh38" | Genome == "hg19"){
      gene.markers <- toupper(gene.markers)
    }
    else if (Genome == "mm10"){
      gene.markers <- str_to_title(gene.markers)
    }
    else if (Genome == "Dualhg19"){
      gene.markers <- paste("hg19-", toupper(gene.markers), sep="")
    }
    else if (Genome == "DualGRCh38"){
      gene.markers <- paste("GRCh38-", toupper(gene.markers), sep="")
    }
    else if (Genome == "Dualmm10"){
      gene.markers <- paste("mm10---", str_to_title(gene.markers), sep="")
    }
    
    gene.markers <- gene.markers[gene.markers %in% getGenes(Seurat.obj)]
    
    lapply(X = gene.markers, FUN = function(x)
    {
      print(FeaturePlot(Seurat.obj, features = x, label = TRUE))
    })
    
  }
  
  print(DimPlot(Seurat.obj, group.by = 'Phase'))
  print(DimPlot(Seurat.obj, group.by = "ID"))
  dev.off()
  
  return(Seurat.obj)
}
############################################################################################################





############################################################################################################
Calculate.Cell.Type.Signature <- function(Seurat.obj, Project.Path, Figure.Total, Analysis.Name, Integration.Method, Cell.Type.DF, Genome)
{
  pdf(file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_all_cell_type_signatures.pdf", sep = ""), width = length(Cell.Type.DF)*6, height = 6)
  
  for (i in 1:ncol(Cell.Type.DF))
  {
    cell.type <- paste(colnames(Cell.Type.DF)[i], ".score", sep = "")
    gene.markers <- as.character(Cell.Type.DF[, i][Cell.Type.DF[, i] != ""])
    
    if (Genome == "GRCh38" | Genome == "hg19"){
      gene.markers <- toupper(gene.markers)
    }
    else if (Genome == "mm10"){
      gene.markers <- str_to_title(gene.markers)
    }
    else if (Genome == "Dualhg19"){
      gene.markers <- paste("hg19-", toupper(gene.markers), sep="")
    }
    else if (Genome == "DualGRCh38"){
      gene.markers <- paste("GRCh38-", toupper(gene.markers), sep="")
    }
    else if (Genome == "Dualmm10"){
      gene.markers <- paste("mm10---", str_to_title(gene.markers), sep="")
    }
    
    gene.markers <- gene.markers[gene.markers %in% getGenes(Seurat.obj)]
    
    Seurat.obj <- AddModuleScore(Seurat.obj, features = list(gene.markers), name = cell.type)
    
    print(FeaturePlot(Seurat.obj, features = paste(cell.type, "1", sep=""), label = TRUE) + theme(aspect.ratio = 1))
    print(VlnPlot(Seurat.obj, features = paste(cell.type, "1", sep=""), pt.size = 0))
    print(VlnPlot(Seurat.obj, features = paste(cell.type, "1", sep=""), pt.size = 0, group.by = "ID"))
  }
  
  
  cell_cluster_scores <- Seurat.obj@meta.data[, grepl("score", colnames(Seurat.obj@meta.data))]
  
  # Predict.Cell.Score <- function(Data.Frame)
  # {
  #   Data.Frame <- sort(Data.Frame, decreasing = TRUE)
  #   if (foldchange(Data.Frame[1], Data.Frame[2]) > 1)
  #   {
  #     celltype <- sub(".score1", "", names(Data.Frame[1]))
  #   }
  #   else
  #   {
  #     celltype <- "UNDETERMINED"
  #   }
  #   
  #   return(celltype)
  # }
  # Seurat.obj$predicted.cell.signature.ident <- apply(cell_cluster_scores, MARGIN = 1, FUN = Predict.Cell.Score)
  
  Seurat.obj$predicted.cell.signature.ident <- sub(".score1", "", colnames(cell_cluster_scores)[max.col(cell_cluster_scores, ties.method = "first")])
  
  print(DimPlot(Seurat.obj, reduction = "umap", group.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
  print(DimPlot(Seurat.obj, reduction = "umap", group.by = "Phase") + theme(aspect.ratio = 1))
  print(DimPlot(Seurat.obj, reduction = "umap", group.by = "ID") + theme(aspect.ratio = 1))
  print(DimPlot(Seurat.obj, reduction = "umap", split.by = "predicted.cell.signature.ident") + theme(aspect.ratio = 1))
  print(ggplot(data = data.frame(table(Seurat.obj$predicted.cell.signature.ident)), aes(x = Var1, y = Freq)) + geom_bar(stat = "identity") + geom_text(aes(label = Freq)))
  
  dev.off()
  
  
  write.table(table(Seurat.obj$predicted.cell.signature.ident, Seurat.obj$ID), file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_all_cell_type_signatures_count_data.tsv", sep = ""), sep="\t", quote=FALSE)
  
  write.table(prop.table(table(Seurat.obj$predicted.cell.signature.ident, Seurat.obj$ID), margin = 2) *100, file = paste(Project.Path, "/", Analysis.Name, "/", Analysis.Name, "_", Figure.Total, "_", Integration.Method, "/", Analysis.Name, "_all_cell_type_signatures_percentage_data.tsv", sep = ""), sep="\t", quote=FALSE)
  
  return(Seurat.obj)
}
############################################################################################################