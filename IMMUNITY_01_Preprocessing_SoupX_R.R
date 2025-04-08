# Load required libraries
library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(knitr)
library(Matrix)  # for writeMM

options(future.globals.maxSize = 2048 * 2048^2)  # sets limit to 800 MB (adjust as needed)

# Define the base sample directory and sample list
sample_directory <- "C:/Users/wwspa/10X_Single_cell_Data/Annotation_siinfek_360_mice/GEX_VDJ/GEX/dp_mult"
# Remove trailing slash if present
sample_directory <- sub("/+$", "", sample_directory)

samples <- c("1_dp_mult", "2_dp_mult", "3_dp_mult", "4_dp_mult")

# Define an output directory for the exported adjusted count matrices
export_dir <- file.path(sample_directory, "adjusted_counts")
if (!dir.exists(export_dir)) {
  dir.create(export_dir)
}

# List the immediate subdirectories in sample_directory
all_folders <- list.dirs(sample_directory, full.names = TRUE, recursive = FALSE)

# Loop over each sample
for (sample_name in samples) {
  message("Processing sample: ", sample_name)
  
  # Identify the folder that starts with the sample name (e.g., "1_CD44hi_st")
  matches <- all_folders[ grepl(paste0("^", sample_name, "(_|$)"), basename(all_folders)) ]
  
  if (length(matches) == 0) {
    message("No folder found starting with ", sample_name, " in ", sample_directory)
    next
  }
  
  # Use the first matching folder and then the 'outs' subfolder
  correct_path <- file.path(matches[1], "outs")
  
  if (!dir.exists(correct_path)) {
    message("Outs directory not found for sample: ", sample_name, " at ", correct_path)
    next
  }
  
  # Load the filtered and raw data using Seurat's Read10X
  filtered_data <- Seurat::Read10X(file.path(correct_path, "filtered_feature_bc_matrix"))
  raw_data <- Seurat::Read10X(file.path(correct_path, "raw_feature_bc_matrix"))
  
  # If the data contains multiple assays, select "Gene Expression" if available
  if ("Gene Expression" %in% names(filtered_data)) {
    filtered_data <- filtered_data[["Gene Expression"]]
  }
  if ("Gene Expression" %in% names(raw_data)) {
    raw_data <- raw_data[["Gene Expression"]]
  }
  
  # Create a SoupChannel object using raw and filtered data
  sc <- SoupChannel(raw_data, filtered_data)
  
  # Create a Seurat object from the filtered data for downstream analysis
  srat <- CreateSeuratObject(counts = filtered_data)
  
  # Run the standard Seurat workflow
  srat <- SCTransform(srat, verbose = FALSE)
  srat <- RunPCA(srat, verbose = FALSE)
  srat <- RunUMAP(srat, dims = 1:30, verbose = FALSE)
  srat <- FindNeighbors(srat, dims = 1:30, verbose = FALSE)
  srat <- FindClusters(srat, verbose = TRUE)
  
  # Extract metadata and UMAP embeddings from the Seurat object
  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  
  # Update the SoupChannel with clustering and dimensionality reduction information
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc <- setDR(sc, umap)
  
  # Estimate contamination and adjust the counts using SoupX
  sc <- autoEstCont(sc)
  adj_matrix <- adjustCounts(sc, roundToInt = TRUE)
  
  # Export the adjusted count matrix and its annotations:
  # Define file paths for the export
  matrix_file <- file.path(export_dir, paste0(sample_name, "_adjusted_counts.mtx"))
  genes_file <- file.path(export_dir, paste0(sample_name, "_genes.tsv"))
  barcodes_file <- file.path(export_dir, paste0(sample_name, "_barcodes.tsv"))
  
  # Write the adjusted count matrix in Matrix Market format
  writeMM(adj_matrix, file = matrix_file)
  
  # Write gene names (rownames) and cell barcodes (colnames)
  write.table(rownames(adj_matrix), file = genes_file, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(colnames(adj_matrix), file = barcodes_file, sep = "\t", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  message("Finished processing sample: ", sample_name)
  message("Exported adjusted counts to:")
  message("  Matrix: ", matrix_file)
  message("  Genes: ", genes_file)
  message("  Barcodes: ", barcodes_file)
}

