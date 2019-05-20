library(Seurat)
library(dplyr)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = snakemake@params[[1]])
# pbmc.data is a matrix of size 32,738 (genes) by 2,700 (cells)

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size

sparse.size <- object.size(x = pbmc.data)
sparse.size

dense.size/sparse.size
# definitely want to use their sparse matrices

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,
                           project = "10X_PBMC")
# resulting matrix is 13,714 genes by 2700 cells

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
# what is the percent of reads that come from mitochondria genes from each cell, which is represented by a unique molecular identifier (UMI)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pdf(snakemake@output[["qc_volcano_plot"]])
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
# nUMI represents the number of reads in a cell

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
# now we are down to 2638 samples

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",
                      scale.factor = 10000)

pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# identifies 1838 genes with high, variable expression
length(x = pbmc@var.genes)

# scaling the data and removing unwanted sources of variation
# this is a basic approach that only regresses on the number of detected molecules per cell as well as the percentage mitochondrial gene content
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# Perform linear dimensional reduction
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,
               genes.print = 5)
# different approaches for visualizing PCA
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
# ProjectPCA scores each gene in the dataset (including genes not included
# in the PCA) based on their correlation with the calculated components.
# Though we don't use this further here, it can be used to identify markers
# that are strongly correlated with cellular heterogeneity, but may not have
# passed through variable gene selection.  The results of the projected PCA
# can be explored by setting use.full=T in the functions above
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)

# Determine statistically significant principal components
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
# pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
# JackStrawPlot(object = pbmc, PCs = 1:12)
PCElbowPlot(object = pbmc)

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10,
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)

# Run non-linear dimensional reduction (tSNE)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = pbmc)
saveRDS(pbmc, file = snakemake@output[["seurat_object"]])
