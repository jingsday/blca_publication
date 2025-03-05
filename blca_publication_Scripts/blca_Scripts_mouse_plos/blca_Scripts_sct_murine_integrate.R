# Load the required packages
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)

#Part I: Read data, merge, create Seurat object and preprocessing(if provided raw)

setwd('/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_mouse_GSE174182_RAW')

names_list <- c("GSM5288668_Sample-3_", "GSM5288669_Sample-4_","GSM5288670_Sample-5_" ,
                "GSM5288671_Sample-6_", "GSM5288672_Sample-7_", "GSM5288673_Sample-8_",
                "GSM5288674_Sample-11_")


for (name in names_list) {
  cts <- ReadMtx(mtx = paste0(name,'filtered_matrix.mtx.gz'),
                 features = paste0(name,'filtered_features.tsv.gz'),
                 cells = paste0(name,'filtered_barcodes.tsv.gz'))
  
  # create seurat objects
  assign(str_sub(name, end = 10),CreateSeuratObject(counts = cts))
}

ls()
merged_seurat <- merge(GSM5288668, y = c(GSM5288669,GSM5288670, GSM5288671, GSM5288672,GSM5288673, GSM5288674),
                       add.cell.ids = ls()[2:8],
                       project = 'BLCA')
#before it was Mt
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^mt-')

# create a sample column
merged_seurat$sample <-rownames(merged_seurat@meta.data)

# split sample column
merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Sample', 'Barcode'), 
                                    sep = '_')

#Murine cells were filtered to retain higher quality cells (>200 &<8000 uniquely identified genes
#<25% of reads mapped to mitochondrial genes),
merged_seurat_filtered <- subset(merged_seurat, subset = nFeature_RNA > 200 &
                                   nFeature_RNA < 8000 &
                                   mitoPercent < 25)



obj.list <- SplitObject(merged_seurat_filtered, split.by = 'Sample')

table(merged_seurat_filtered@meta.data$Sample)

#Part II: SCT intergration and save files for furthur analysis

obj.list <- lapply(X = obj.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                  anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

#Save intergrated data 
outdir <- '/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/'
# saveRDS(seurat.integrated, file = paste0(outdir,"murine_intergration_sct.rds"))

seurat.integrated <-readRDS('/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/blca_publication_OUTPUT_sct/murine_intergration_sct.rds')
#Downstream clustering, visualization and find markers
DefaultAssay(seurat.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
seurat.integrated <- RunPCA(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated,  dims = 1:30)

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30, verbose = FALSE)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.05,random.seed =1)

table(seurat.integrated@meta.data$seurat_clusters)#,seurat.integrated@meta.data$Sample)

markers <- c('Cdh1', 'Upk1a', 'Upk1b', 'Upk2', 'Upk3a', 'Ivl','Krt7',"Krt8", "Krt18", "Krt20","Gata3", "Foxa1")
# Urothelial cell markers for mouse bladder
markers <- c('Cdh1', 'Upk1a', 'Upk2', 'Upk3a', 'Ivl','Krt7', "Foxa1")
plot.new()
DotPlot(seurat.integrated, features = markers) + RotatedAxis()

png(paste0(outdir,'blca_publication_OUTPUT_fig/blca_publication_fig_markers.png'),width = 2400,height = 1500,res=300)
DotPlot(seurat.integrated, features = markers) + RotatedAxis()
dev.off()
# Single cell heatmap of feature expression
DoHeatmap(subset(seurat.integrated, downsample = 100), features = markers, size = 3)

# 
# seurat.integrated.markers <- FindAllMarkers(seurat.integrated, only.pos = TRUE)
# seurat.integrated.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1 & p_val_adj <0.05)
# 
# seurat.integrated.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 5) %>%
#   ungroup() -> top10
# 
# png(paste0(outdir,'blca_publication_fig_markers.png'),res=300)
# DoHeatmap(seurat.integrated, features = top10$gene) + NoLegend()
# dev.off()

#visualization
p2 <- DimPlot(seurat.integrated, reduction = "umap")
p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by = "Sample", 
              repel = TRUE,cols=c('red','green','blue','yellow','brown','orange','purple'))
p1 + p2
png(paste0(outdir,'blca_publication_OUTPUT_fig/blca_publication_source.png'),width = 2400,height = 1500,res=300)
p1 + p2
dev.off()

library(ggplot2)
library(dplyr)

# Extract metadata
df <- seurat.integrated@meta.data %>%
  filter(seurat_clusters == 0) %>%  # Only keep cluster 0
  count(Sample) %>%  
  mutate(Frequency = sum(n))  # Compute relative frequency

df$Sample <- factor(df$Sample, levels = c("GSM5288674", "GSM5288673", "GSM5288672",
                                          "GSM5288671", "GSM5288670", "GSM5288669", 
                                          "GSM5288668"))  
# Horizontal Bar Plot (Cell Counts)
ggplot(df, aes(x = Sample, y = n, fill = Sample)) +  # Order by counts
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Add border color
  scale_fill_manual(values = c('red', 'green', 'blue', 'yellow', 'brown', 'orange', 'purple')) +  # Assign custom colors
  coord_flip() +  # Make bars horizontal
  labs(x = "Sample", y = "Cell Counts", fill = "Sample") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14,),  # Adjust y-axis text size
    axis.text.x = element_text(size = 14,angle = 90),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
png(paste0(outdir,'blca_publication_OUTPUT_fig/blca_publication_fig_bar.png'),width = 2200,height = 1800,res=300)
dev.off()


ggplot(df, aes(x = "", y = Frequency, fill = Sample)) +
  scale_fill_manual(values = c('red', 'green', 'blue', 'yellow', 'brown', 'orange', 'purple')) +  # Assign custom colors
  
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add border color if needed
  labs(fill = "Sample") +
  theme_void() +  
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )

# Define Sample order
df$Sample <- factor(df$Sample, levels = c("GSM5288674", "GSM5288673", "GSM5288672",
                                          "GSM5288671", "GSM5288670", "GSM5288669", 
                                          "GSM5288668"))  

# Plot pie chart
ggplot(df, aes(x = "", y = Frequency, fill = Sample)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # Add border color if needed
  coord_polar(theta = "y") +  # Convert to pie chart
  scale_fill_manual(values = c('red', 'green', 'blue', 'yellow', 'brown', 'orange', 'purple')) +  # Assign custom colors
  labs(fill = "Sample") +
  theme_void() +  
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
png(paste0(outdir,'blca_publication_OUTPUT_fig/blca_publication_fig_pie.png'),width = 2200,height = 1800,res=300)
dev.off()


df <- seurat.integrated@meta.data %>%
  count(Sample, seurat_clusters==0) %>%  # Count cells per cluster per sample
  group_by(Sample) %>%
  mutate(Frequency = n / sum(n))  # Compute relative frequency
df$Sample <- factor(df$Sample, levels = c("GSM5288674", "GSM5288673", "GSM5288672" ,"GSM5288671", "GSM5288670", "GSM5288669" ,"GSM5288668"))  # Specify your order
p <- ggplot(df, aes(x = Sample, y = Frequency, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # Make bars horizontal
  labs(y='Cell Type Fraction', fill = "Louvain Cluster") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),  # Align text sizes
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
    )
  )
#Double check
table(seurat.integrated@meta.data$seurat_clusters,seurat.integrated@meta.data$Sample)
p
#Save
png(paste0(outdir,'blca_publication_OUTPUT_fig/blca_publication_fig_fraction.png'),width = 2000,height = 1800,res=300)
p
dev.off()


#Visualize marker genes as violin plots.
DefaultAssay(seurat.integrated) <- "RNA"

p1 <- VlnPlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"),
              pt.size = 0.2, ncol = 4)
p1

p1 <- FeaturePlot(seurat.integrated, features = c("Cdh1", "Ivl", "Upk1a", "Upk1b", "Upk3a", "Upk3b", "Slc14a1", "Upk2"), pt.size = 0.2,
                  ncol = 3)
p1+p2

#Retrieve sct data 
sct_murine_cl00 <- subset(seurat.integrated, subset = seurat_clusters == 0)

corrected_UMI <- sct_murine_cl00[["SCT"]]$data
corrected_UMI <-t(corrected_UMI)
writeMM(corrected_UMI, paste0(outdir,"blca_publication_OUTPUT_sct/sct_corrected_UMI.mtx"))
write.table(rownames(corrected_UMI), paste0(outdir, "blca_publication_OUTPUT_sct/sct_corrected_UMI_cells.txt"), col.names = FALSE, row.names = FALSE)
write.table(colnames(corrected_UMI), paste0(outdir,"blca_publication_OUTPUT_sct/sct_corrected_UMI_genes.txt"), col.names = FALSE, row.names = FALSE)



colnames_substr <- substr(colnames(corrected_UMI), 1, 10)

#input interested datasets
indices <- which(colnames_substr == "GSM5288674")
GSM5288674_corrected_UMI <- corrected_UMI[, indices]
GSM5288674_corrected_UMI<- t(GSM5288674_corrected_UMI)
writeMM(GSM5288674_corrected_UMI, "sct_full_UM/GSM5288674_corrected_UMI.mtx")
write.table(rownames(GSM5288674_corrected_UMI), file = "sct_full_UM/GSM5288674_corrected_UMI_cells.txt", col.names = FALSE, row.names = FALSE)
write.table(colnames(GSM5288674_corrected_UMI), file = "sct_full_UM/GSM5288674_corrected_UMI_genes.txt", col.names = FALSE, row.names = FALSE)

