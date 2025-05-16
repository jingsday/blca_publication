library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)
data_dir <- '/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_GSE135337/'
dirs <- list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz

names_list <- c()  # Create an empty list to store the names
for (x in dirs) {
  name <- unique(str_sub(x, end = 14))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}
names_list <-unique(names_list)
names_list
#Select Ta low grade as they don't invade 
names_list <- names_list[c(1,4,5,6,7)]
names_list

for (name in names_list) {
  cts <- read.delim(paste0(data_dir, name, "_gene_cell_exprs_table.txt.gz"))
  cts[, 2]<- make.unique(cts[, 2])
  
  rownames(cts) <- cts[, 2]  #first col ensemble ID, 2nd gene name
  cts <- cts[, -c(1, 2)]
  
  colnames(cts) <- colnames(cts)  
  assign(name, CreateSeuratObject(counts = cts,min.features = 200))
}

ls()
#1,6,7  Ta  4,5 T2 and T3  
merged_seurat <- merge(GSM4006644_BC1, y = c(GSM4006647_BC4,GSM4006648_BC5, GSM4751267_BC6, GSM4751268_BC7),
                       add.cell.ids = ls()[4:8],
                       project = 'BLCA')
#Ta vs MIBC  but need to read. Luminal T3, even mixed cell types in Ta
str(merged_seurat)


merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample','Patient', 'Barcode'), 
                                    sep = '_')
head(merged_seurat@meta.data)


#Checking QC 
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


merged_seurat_qc <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6001)
table(merged_seurat_qc@meta.data$Patient)
dim(merged_seurat_qc@meta.data)
table(merged_seurat@meta.data$Sample)
dim(merged_seurat@meta.data)

#scTransform 
obj.list <- SplitObject(merged_seurat_qc, split.by = 'Sample')

obj.list <- lapply(X = obj.list, FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)

obj.list <- PrepSCTIntegration(
  obj.list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

outdir <- '/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/'
#saveRDS(seurat.integrated, file = paste0(outdir,"human_intergration_sct.rds"))
#worked this time 

corrected_UMI <- seurat.integrated[["SCT"]]$data
corrected_UMI <-t(corrected_UMI)

writeMM(corrected_UMI, paste0(outdir,"blca_scR_corrected_UMI.mtx"))
write.table(rownames(corrected_UMI), file = paste0(outdir, "blca_scR_corrected_UMI_cells.txt"), col.names = FALSE, row.names = FALSE)
write.table(colnames(corrected_UMI), file = paste0(outdir, "blca_scR_corrected_UMI_genes.txt"), col.names = FALSE, row.names = FALSE)

#check gene row length. Might be mismatch. Version 1 or v2 might solve. 
#One more thing 

###DEG 

seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
seurat.integrated <- FindNeighbors(seurat.integrated, reduction = "pca", dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.3)


DimPlot(seurat.integrated, reduction = "umap")
head(seurat.integrated@meta.data)
#Passing annnotation 
annotation <- read.delim('/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/non_immune_cells_obs.csv',sep=',')
annotation$IDs <- paste0(annotation$source,'_',  annotation$CellID)

seurat_meta <- seurat.integrated@meta.data
seurat_meta$cellID <- rownames(seurat_meta)

# Then merge
merged_df <- merge(seurat_meta, annotation[, c('IDs', 'leiden', 'Stage')],
                   by.x = "cellID", by.y = "IDs", all.x = TRUE)

rownames(merged_df) <- merged_df$cellID
seurat.integrated@meta.data <- merged_df[, setdiff(colnames(merged_df), "cellID")]
head(seurat.integrated@meta.data)

seurat.integrated@meta.data$target_cells <- paste0(seurat.integrated@meta.data$leiden,'_',seurat.integrated@meta.data$Stage)
DimPlot(seurat.integrated, reduction = "umap", split.by = "Stage")

table(seurat.integrated@meta.data$target_cells)

Idents(seurat.integrated) <- "target_cells"
seurat.integrated <- PrepSCTFindMarkers(seurat.integrated)

invasive_sig <- FindMarkers(seurat.integrated, assay = "SCT", ident.1 = "1_T2", ident.2 = "1_Ta",
                                     verbose = FALSE)
head(invasive_sig, n = 15)

rownames(invasive_sig[invasive_sig$p_val_adj < 0.05, ])#6669

#write.csv(invasive_sig,paste0(outdir,'invasive_sig.csv'))
#Check overlap with LINCS genes 
library(readxl)

stv <- read_excel('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_LINCS/ALL_DATA_2020_Jing_gbm_del.xlsx',
                  sheet = 'STVs')
common_lincs <- intersect(rownames(invasive_sig[invasive_sig$p_val_adj < 0.05, ]), stv$Gene) # 532


#
DefaultAssay(seurat.integrated) <-'SCT'
FeaturePlot(seurat.integrated, features = c("CD36", "S100A8", "SPINK1"), split.by = "Stage", max.cutoff = 3,
            cols = c("grey", "red"))
Idents(seurat.integrated) <- "leiden"

cluster1 <- FindConservedMarkers(seurat.integrated, assay = "SCT", ident.1 = "1", grouping.var = "Stage",
                                   verbose = FALSE)
head(cluster1)
FeaturePlot(seurat.integrated, features = c("RPS27", "FTH1", "CSTB"), split.by = "Stage", max.cutoff = 3,
            cols = c("grey", "red"))
VlnPlot(seurat.integrated, features = c("S100A7", "S100A8", "S100A9"), split.by = "Stage",
                 group.by = "leiden", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
cluster1 <- cluster1[cluster1$max_pval<0.05,]
#write.csv(cluster1,paste0(outdir,'canonical_c1_sig.csv'))

