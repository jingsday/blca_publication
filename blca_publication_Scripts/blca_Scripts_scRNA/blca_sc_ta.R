library(Seurat)

# For output from CellRanger >= 3.0 with multiple data types
data_dir <- '/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_GSE135337/'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir,gene.column = 1,)


seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
seurat_object[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)


counts <- read.delim(paste0(data_dir,"GSM5329919_BCN_gene_cell_exprs_table.xls.gz"), row.names=1)

str(counts)

seurat_object = CreateSeuratObject(counts = counts)

str(seurat_object)

data <- read.delim("~/Desktop/GSE155512_RAW/GSM4705589_RPE004_matrix.txt.gz", header = T, stringsAsFactors = F)
seurat <- CreateSeuratObject(data, min.cells = 10, min.features = 200, project = "seurat")
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

#Ta vs MIBC  but need to read. Luminal T3, even mixed cell types in Ta