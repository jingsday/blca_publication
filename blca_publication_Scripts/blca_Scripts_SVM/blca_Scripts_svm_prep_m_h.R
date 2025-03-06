library(biomaRt)
library(dplyr)

datadir <-'/home/jing/Phd_project/project_UCD_blca/blca_publication_OUTPUT/'

# Load mouse gene names
features <- read.csv(paste0(datadir, 'blca_publication_OUTPUT_sct/', 'sct_corrected_UMI_genes.txt'), 
                     header = FALSE, stringsAsFactors = FALSE)

# Create a dataframe with gene names
name_df <- data.frame(mouse_gene = features$V1)  # Assuming gene names are in the first column
name_df$orig.id <- name_df$mouse_gene

# Convert to character type (to prevent errors in joins)
name_df$mouse_gene <- as.character(name_df$mouse_gene)

# Connect to Ensembl and select the mouse dataset
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Query Ensembl to get human homologs
gene_IDs <- getBM(
  filters = "external_gene_name",
  attributes = c("external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  values = name_df$mouse_gene,
  mart = mouse_mart
)

# Merge mouse genes with human homologs
gene.df <- left_join(name_df, gene_IDs, by = c("mouse_gene" = "external_gene_name"))

# Remove rows without human homologs
gene.df <- gene.df[!is.na(gene.df$hsapiens_homolog_associated_gene_name), ]

# Ensure row names are unique
rownames(gene.df) <- make.unique(gene.df$orig.id)


write.csv(gene.df,paste0(datadir, 'blca_publication_OUTPUT_sct/','blca_publication_OUTPUT_m_h.csv'))

          