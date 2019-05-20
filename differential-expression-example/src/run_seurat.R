library(Seurat)

# Functions
LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0,
    adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
    for (i in genes) {
        x1 <- exp.mat[i, 1]
        y1 <- exp.mat[i, 2]
        plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t,
            label = i, size = text.size)
        plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 +
            adj.y.s, yend = y1, size = segment.size)
    }
    return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05,
    adj.r.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t,
        adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05,
    adj.l.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t,
        adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}

# main
ctrl.data <- read.table(snakemake@input[["ctrl"]], sep = "\t")
stim.data <- read.table(snakemake@input[["stim"]], sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(raw.data = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)
# Set up stimulated object
stim <- CreateSeuratObject(raw.data = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim@meta.data$stim <- "STIM"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)

# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))

# Perform canonical correlation analysis
immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = immune.combined, reduction.use = "cca", group.by = "stim",
    pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "CC1", group.by = "stim",
    do.return = TRUE)
plot_grid(p1, p2)

# use this statistic to estimate how many correlation components (CC) to use for downstream analysis
p3 <- MetageneBicorPlot(immune.combined, grouping.var = "stim", dims.eval = 1:30,
                        display.progress = FALSE)

# visualize the genes driving each CC
DimHeatmap(object = immune.combined, reduction.type = "cca", cells.use = 500,
           dim.use = 1:9, do.balanced = TRUE)

# align the CCA subspaces - why? batch correction?
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim",
    dims.align = 1:20)

p1 <- VlnPlot(object = immune.combined, features.plot = "ACC1", group.by = "stim",
    do.return = TRUE)
p2 <- VlnPlot(object = immune.combined, features.plot = "ACC2", group.by = "stim",
    do.return = TRUE)
plot_grid(p1, p2)

# t-SNE and Clustering
immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20,
    do.fast = T)
immune.combined <- FindClusters(immune.combined, reduction.type = "cca.aligned",
    resolution = 0.6, dims.use = 1:20)

# Visualization
p1 <- TSNEPlot(immune.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(immune.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

# Identify cell type markers that are conserved across both conditions
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim",
    print.bar = FALSE) # specifically for cluster 7, which happens to be NK cells
head(nk.markers)

# manually pick out conserved cell type markers and use them to visualize the clusters
FeaturePlot(object = immune.combined, features.plot = c("CD3D", "SELL", "CREM",
    "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9", cols.use = c("lightgrey",
    "blue"), pt.size = 0.5)

# label the clusters using our manually identified labels
new.ident <- c("CD14 Mono", "CD4 Naive T", "CD4 Memory T", "B", "CD16 Mono",
    "T activated", "CD8 T", "NK", "DC", "B activated", "Mk", "pDC", "Eryth")
for (i in 0:12) {
    immune.combined <- RenameIdent(object = immune.combined, old.ident.name = i,
        new.ident.name = new.ident[i + 1])
}

TSNEPlot(immune.combined, do.label = T, pt.size = 0.5)

# SplitDotPlotGG is useful for viewing conserved cell type markers across conditions
immune.combined@ident <- factor(immune.combined@ident, levels = (c("pDC", "Eryth",
    "Mk", "DC", "CD14 Mono", "CD16 Mono", "B activated", "B", "CD8 T", "NK",
    "T activated", "CD4 Naive T", "CD4 Memory T")))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY",
    "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A",
    "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2",
    "HBB", "TSPAN13", "IL3RA", "IGJ")
sdp <- SplitDotPlotGG(immune.combined, genes.plot = rev(markers.to.plot), cols.use = c("blue",
    "red"), x.lab.rot = T, plot.legend = T, dot.scale = 8, do.return = T, grouping.var = "stim")

# Identify genes that are differentially expressed across conditions (in the same cell type)

# Focus on naive T cells and CD14 monocyte populations and plot the differences
t.cells <- SubsetData(immune.combined, ident.use = "CD4 Naive T", subset.raw = T)
t.cells <- SetAllIdent(t.cells, id = "stim")
avg.t.cells <- log1p(AverageExpression(t.cells, show.progress = FALSE))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- SubsetData(immune.combined, ident.use = "CD14 Mono", subset.raw = T)
cd14.mono <- SetAllIdent(cd14.mono, id = "stim")
avg.cd14.mono <- log1p(AverageExpression(cd14.mono, show.progress = FALSE))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label1 = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1")
genes.to.label2 = c("IFIT2", "IFIT1")
genes.to.label3 = c("CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avg.t.cells,
    adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avg.t.cells, adj.u.t = 0.5, adj.u.s = 0.4,
    adj.l.t = 0.25, adj.l.s = 0.25)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avg.cd14.mono,
    adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avg.cd14.mono, adj.u.t = 0.5, adj.u.s = 0.4,
    adj.l.t = 0.25, adj.l.s = 0.25)
plot_grid(p1, p2)

# now focus on B cells and print the statistics
immune.combined@meta.data$celltype.stim <- paste0(immune.combined@ident, "_",
    immune.combined@meta.data$stim)
immune.combined <- StashIdent(immune.combined, save.name = "celltype")
immune.combined <- SetAllIdent(immune.combined, id = "celltype.stim")
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL",
    print.bar = FALSE)
head(b.interferon.response, 15)

# good visualization for the change in specific genes across all clusters in the two conditions
FeatureHeatmap(immune.combined, features.plot = c("CD3D", "GNLY", "IFI6", "ISG15",
    "CD14", "CXCL10"), group.by = "stim", pt.size = 0.25, key.position = "top",
    max.exp = 3)

saveRDS(immune.combined, file = snakemake@output[["seurat_file"]])
