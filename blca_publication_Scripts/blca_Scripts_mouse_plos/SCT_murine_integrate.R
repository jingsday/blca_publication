# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)

#Part I: Read data, merge, create Seurat object and preprocessing(if provided raw)

setwd('/home/jing/resource/murine_blca')
outdir <- '/home/jing/resource/murine_blca/scripts/'

# get data location
dirs <- list.files(path = '/home/jing/resource/murine_blca/GSE174182_RAW', recursive = F, full.names = F)


names_list <- c()  # Create an empty list to store the names
for (x in dirs) {
  name <- unique(str_sub(x, end = 20))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}
names_list <-unique(names_list)

names_list <- names_list[names_list != "GSM5288674_Sample-11"]

# Add "GSM5288674_Sample-11_" to the list
names_list <- c(names_list, "GSM5288674_Sample-11_")


for (name in names_list) {
  cts <- ReadMtx(mtx = paste0('GSE174182_RAW/',name,'filtered_matrix.mtx.gz'),
                 features = paste0('GSE174182_RAW/',name,'filtered_features.tsv.gz'),
                 cells = paste0('GSE174182_RAW/',name,'filtered_barcodes.tsv.gz'))
  
  # create seurat objects
  assign(str_sub(name, end = 10), CreateSeuratObject(counts = cts))
}



ls()

merged_seurat <- merge(GSM5288668, y = c(GSM5288669,GSM5288670, GSM5288671, GSM5288672, GSM5288673, GSM5288674),
                       add.cell.ids = ls()[3:9],
                       project = 'BLCA')

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

obj.list <- SplitObject(merged_seurat, split.by = 'Sample')

#Part II: SCT intergration and save files for furthur analysis

obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#Save intergrated data 
#saveRDS(seurat.integrated, file = "murine_intergration_sct.rds")

#Barcodes and samples of interest
cluster_0_cells <- which(seurat.integrated@meta.data$seurat_clusters == 0)

# Subset seurat.integrated to include only cells belonging to cluster 0
seurat_integrated_cluster_0 <- seurat.integrated[cluster_0_cells, ]


sct_murine_cl00 <- subset(sct_murine, subset = seurat_clusters == 0)
saveRDS(sct_murine_cl00, file = "sct_murine_cl00.rds")
#sct_murine_cl00[['SCT]] for processing ultimately not here

#Downstream clustering, visualization and find markers
DefaultAssay(seurat.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.5)


head(seurat.integrated@meta.data@seurat)

#visualization

p1 <- DimPlot(seurat.integrated, reduction = "umap")
p2 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Sample", 
              repel = TRUE,cols=c('red','green','blue','yellow','brown','orange','purple'))

p1 + p2
combined <- p1+p2

#Visualize marker genes as violin plots.
DefaultAssay(seurat.integrated) <- "RNA"

p1 <- VlnPlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"),
              pt.size = 0.2, ncol = 4)

p1 <- FeaturePlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"), pt.size = 0.2,
                  ncol = 3)
p1+p2

