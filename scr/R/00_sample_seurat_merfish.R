# Analysis of Image-based Spatial Data in Seurat
# Compiled: May 02, 2022
# Source: vignettes/spatial_vignette_2.Rmd
# Overview
# In this vignette, we introduce a Seurat extension to analyze new types of spatially-resolved data. We have previously introduced a spatial framework which is compatible with sequencing-based technologies, like the 10x Genomics Visium system, or SLIDE-seq. Here, we extend this framework to analyze new data types that are captured via highly multiplexed imaging. In contrast to sequencing-based technologies, these datasets are often targeted (i.e. they profile a pre-selected set of genes). However they can resolve individual molecules - retaining single-cell (and subcellular) resolution. These approaches also often capture cellular boundaries (segmentations).
# 
# We update the Seurat infrastructure to enable the analysis, visualization, and exploration of these exciting datasets. In this vignette, we focus on three datasets produced by different multiplexed imaging technologies, each of which is publicly available. We will be adding support for additional imaging-based technologies in the coming months.
# 
# Vizgen MERSCOPE (Mouse Brain)
# Nanostring CosMx Spatial Molecular Imager (FFPE Human Lung)
# Akoya CODEX (Human Lymph Node)
# First, we install the updated versions of Seurat and SeuratObject that support this infrastructure, as well as other packages necessary for this vignette.

# library(remotes)
# remotes::install_github("satijalab/seurat", "feat/imaging")
library(Seurat)
library(future)
library(hdf5r)
library(tidyverse)
plan("multisession", workers = 10)

# Mouse Brain: Vizgen MERSCOPE --------------------------------------------
# This dataset was produced using the Vizgen MERSCOPE system, which utilizes the MERFISH technology. The total dataset is available for public download, and contains nine samples (three full coronal slices of the mouse brain, with three biological replicates per slice). The gene panel consists of 483 gene targets, representing known anonical cell type markers, nonsensory G-Protein coupled receptors (GPCRs), and Receptor Tyrosine Kinases (RTKs). In this vignette, we analyze one of the samples - slice 2, replicate 1. The median number of transcripts detected in each cell is 206.

# First, we read in the dataset and create a Seurat object.

# We use the LoadVizgen() function, which we have written to read in the output of the Vizgen analysis pipeline. The resulting Seurat object contains the following information:
#   
#   A count matrix, indicating the number of observed molecules for each of the 483 transcripts in each cell. This matrix is analogous to a count matrix in scRNA-seq, and is stored by default in the RNA assay of the Seurat object

# Loading segmentations is a slow process and multi processing with the future pacakge is
# recommended

vizgen.obj <- LoadVizgen(data.dir = "../data/mouse_brain_map/BrainReceptorShowcase/Slice2/Replicate1/", fov = "s2r1")

# The next pieces of information are specific to imaging assays, and is stored in the images slot of the resulting Seurat object:
  
# Cell Centroids: The spatial coordinates marking the centroid for each cell being profiled
# Get the center position of each centroid. There is one row per cell in this dataframe.
head(GetTissueCoordinates(vizgen.obj[["s2r1"]][["centroids"]]))

# Cell Segmentation Boundaries: The spatial coordinates that describe the polygon segmentation of each single cell
# Get the coordinates for each segmentation vertice. Each cell will have a variable number of
# vertices describing its shape.
head(GetTissueCoordinates(vizgen.obj[["s2r1"]][["segmentation"]]))

# Molecule positions: The spatial coordinates for each individual molecule that was detected during the multiplexed smFISH experiment.
# Fetch molecules positions for Chrm1
head(FetchData(vizgen.obj[["s2r1"]][["molecules"]], vars = "Chrm1"))

# Preprocessing and unsupervised analysis ---------------------------------
# We start by performing a standard unsupervised clustering analysis, essentially first treating the dataset as an scRNA-seq experiment. We use SCTransform-based normalization, though we slightly modify the default clipping parameters to mitigate the effect of outliers that we occasionally observe in smFISH experiments. After normalization, we can run dimensional reduction and clustering.
vizgen.obj <- SCTransform(vizgen.obj, assay = "Vizgen", clip.range = c(-10, 10), ) %>%
  RunPCA(npcs = 30, features = rownames(vizgen.obj)) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.3)

# print the object
vizgen.obj

# We can then visualize the results of the clustering either in UMAP space (with DimPlot()) or overlaid on the image with ImageDimPlot().
DimPlot(vizgen.obj, reduction = "umap",label = T,repel = T)
ggsave("out/image/UMAP_s2r1.png",width = 8,height = 8)
ggsave("out/image/UMAP_s2r1.pdf",width = 8,height = 8)

ImageDimPlot(vizgen.obj, fov = "s2r1", cols = "polychrome", axes = TRUE)
ggsave("out/image/OverlayCluster_s2r1.png",width = 8,height = 8)
ggsave("out/image/OverlayCluster_s2r1.pdf",width = 8,height = 8)

# You can also customize multiple aspect of the plot, including the color scheme, cell border widths, and size (see below).

# Customizing spatial plots in Seurat
# The ImageDimPlot() and ImageFeaturePlot() functions have a few parameters which you can customize individual visualizations. These include:
#   
# alpha: Ranges from 0 to 1. Sets the transparency of within-cell coloring.
# size: determines the size of points representing cells, if centroids are being plotted
# cols: Sets the color scheme for the internal shading of each cell. Examples settings are polychrome, glasbey, Paired, Set3, and parade. Default is the ggplot2 color palette
# shuffle.cols: In some cases the selection of cols is more effective when the same colors are assigned to different clusters. Set shuffle.cols = TRUE to randomly shuffle the colors in the palette.
# border.size: Sets the width of the cell segmentation borders. By default, segmentations are plotted with a border size of 0.3 and centroids are plotted without border.
# border.color: Sets the color of the cell segmentation borders
# dark.background: Sets a black background color (TRUE by default)
# axes: Display
# Since it can be difficult to visualize the spatial localization patterns of an individual cluster when viewing them all together, we can highlight all cells that belong to a particular cluster:
p1 <- ImageDimPlot(vizgen.obj, fov = "s2r1", cols = "red", cells = WhichCells(vizgen.obj, idents = 14))
p2 <- ImageDimPlot(vizgen.obj, fov = "s2r1", cols = "red", cells = WhichCells(vizgen.obj, idents = 15))
(p1 + p2) & theme(plot.background = element_rect(fill = "black"))
ggsave("out/image/Cluster14-15_s2r1.png",width = 16,height = 8)
ggsave("out/image/Cluster14-15_s2r1.pdf",width = 16,height = 8)

# We can find markers of individual clusters and visualize their spatial expression pattern. We can color cells based on their quantified expression of an individual gene, using ImageFeaturePlot(), which is analagous to the FeaturePlot() function for visualizing expression on a 2D embedding. Since MERFISH images individual molecules, we can also visualize the location of individual molecules.
p12 <- ImageFeaturePlot(vizgen.obj, features = "Slc17a7")
p22 <- ImageDimPlot(vizgen.obj, molecules = "Slc17a7", nmols = 10000, alpha = 0.3, mols.cols = "red")
(p12 + p22) & theme(plot.background = element_rect(fill = "black"))
ggsave("out/image/Overlay-Slc17a7_s2r1.png",width = 16,height = 8)
ggsave("out/image/Overlay-Slc17a7_s2r1.pdf",width = 16,height = 8)

# Note that the nmols parameter can be used to reduce the total number of molecules shown to reduce overplotting. You can also use the mols.size, mols.cols, and mols.alpha parameter to further optimize.

# Plotting molecules is especially useful for visualizing co-expression of multiple genes on the same plot.
p13 <- ImageDimPlot(vizgen.obj, fov = "s2r1", alpha = 0.3, molecules = c("Slc17a7", "Olig1"), nmols = 10000)

# define the top markers of cluster 14 (14 vs all)
markers.14 <- FindMarkers(vizgen.obj, ident.1 = "14")

# plot the top 4 markers of cluster 14 (14 vs all)
p23 <- ImageDimPlot(vizgen.obj, fov = "s2r1", alpha = 0.3, molecules = rownames(markers.14)[1:4],
                   nmols = 10000)
(p13 + p23) & theme(plot.background = element_rect(fill = "black"))
ggsave("out/image/Overlay-GOI_s2r1.png",width = 16,height = 8)
ggsave("out/image/Overlay-GOI_s2r1.pdf",width = 16,height = 8)

# The updated Seurat spatial framework has the option to treat cells as individual points, or also to visualize cell boundaries (segmentations). By default, Seurat ignores cell segmentations and treats each cell as a point (‘centroids’). This speeds up plotting, especially when looking at large areas, where cell boundaries are too small to visualize.

# We can zoom into a region of tissue, creating a new field of view. For example, we can zoom into a region that contains the hippocampus. Once zoomed-in, we can set DefaultBoundary() to show cell segmentations. You can also ‘simplify’ the cell segmentations, reducing the number of edges in each polygon to speed up plotting.

# create a Crop
cropped.coords <- Crop(vizgen.obj[["s2r1"]], x = c(1750, 3000), y = c(3750, 5250), coords = "plot")

# set a new field of view (fov)
vizgen.obj[["hippo"]] <- cropped.coords

# visualize FOV using default settings (no cell boundaries)
p14 <- ImageDimPlot(vizgen.obj, fov = "hippo", axes = TRUE, size = 0.7, border.color = "white", cols = "polychrome",
                   coord.fixed = FALSE)

# visualize FOV with full cell segmentations
DefaultBoundary(vizgen.obj[["hippo"]]) <- "segmentation"

p24 <- ImageDimPlot(vizgen.obj, fov = "hippo", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)

# simplify cell segmentations
vizgen.obj[["hippo"]][["simplified.segmentations"]] <- Simplify(coords = vizgen.obj[["hippo"]][["segmentation"]],
                                                                tol = 3)
DefaultBoundary(vizgen.obj[["hippo"]]) <- "simplified.segmentations"

# visualize FOV with simplified cell segmentations
DefaultBoundary(vizgen.obj[["hippo"]]) <- "simplified.segmentations"
p34 <- ImageDimPlot(vizgen.obj, fov = "hippo", axes = TRUE, border.color = "white", border.size = 0.1,
                   cols = "polychrome", coord.fixed = FALSE)

(p14 + p24 + p34) & theme(plot.background = element_rect(fill = "black"))
ggsave("out/image/Overlay-cropped_s2r1.png",width = 24,height = 8)
ggsave("out/image/Overlay-cropped_s2r1.pdf",width = 24,height = 8)

# What is the tol parameter?
# We can visualize individual molecules plotted at higher resolution after zooming-in

# Since there is nothing behind the segmentations, alpha will slightly mute colors
ImageDimPlot(vizgen.obj, fov = "hippo", molecules = rownames(markers.14)[1:4], cols = "polychrome",
             mols.size = 1, alpha = 0.5, mols.cols = c("red", "blue", "yellow", "green"))
ggsave("out/image/Overlay-cropped2_s2r1.png",width = 8,height = 8)
ggsave("out/image/Overlay-cropped2_s2r1.pdf",width = 8,height = 8)