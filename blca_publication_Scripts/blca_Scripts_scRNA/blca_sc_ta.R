library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(Matrix)
library(stringr)
library(sctransform)
# For output from CellRanger >= 3.0 with multiple data types
data_dir <- '/home/jing/Phd_project/project_UCD_blca/blca_DATA/blca_DATA_GSE135337/'
dirs <- list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz



names_list <- c()  # Create an empty list to store the names
for (x in dirs) {
  name <- unique(str_sub(x, end = 14))  # Extract the name
  names_list <- c(names_list, name)     # Append the name to the list
}
names_list <-unique(names_list)
names_list

names_list <- names_list[c(1,4,5,6,7)]

#Select Ta low grade as they don't invade 


for (name in names_list) {
  cts <- read.delim(paste0(data_dir,name, "_gene_cell_exprs_table.txt.gz"), row.names=1)  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

ls()
#1,6,7  Ta  4,5 T2 and T3  


merged_seurat <- merge(GSM4006644_BC1, y = c(GSM4006647_BC4,GSM4006648_BC5, GSM4751267_BC6, GSM4751268_BC7),
                       add.cell.ids = ls()[3:9],
                       project = 'BLCA')



#Ta vs MIBC  but need to read. Luminal T3, even mixed cell types in Ta