# 10X Analysis Pipeline: Convert Seurat object to matrix
# Author: Koki Sasagawa

library(argparser, quietly = TRUE)
library(Seurat, quietly = TRUE)
options(future.globals.maxSize = 100000000000)   # Override default memory limit on R

seurat_obj_to_matrix_txt = function(seurat_object, assay_type, slot_type, metadata_cols){
  cellID = colnames(seurat_object)  # Get cell ID's 
  metadata = list()
  for (i in metadata_cols){
    metadata[[i]] = as.vector(seurat_object@meta.data[[i]])
  }
  columnID = c()
  # For each cell, create a new columnID with all associated metadata 
  # This will be saved as a tuple so it can be parsed later in python
  for (i in 1:length(cellID)){
    label = c("('CellID: ", cellID[[i]], "'", ", ")
    for (j in metadata_cols){
      label = c(label, "'", j, ": ", metadata[[j]][i], "'", ", ")
    }
    label = label[-length(label)]  # Remove the last ',' 
    label = c(label, ")")
    label = paste(label, collapse = '')
    columnID = c(columnID, label)
  }
  df = as.data.frame(as.matrix(GetAssayData(seurat_object, slot = slot_type, assay = assay_type)))
  colnames(df) = columnID
  return(df)
}

# Parse command line arguments
p = arg_parser("Convert Seurat object to matrix where rows are genes and columns 
               are cell ID (with optional metadata)")
p = add_argument(p, "--input", help = "Path to seurat object")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--fname", 
                 help = "Name of output file (optional: default is input file name)",
                 default = NA)
p = add_argument(p, "--assay", help = "Assay to access from seurat object")
p = add_argument(p, "--slot", help = "Slot to access from seurat object assay")
p = add_argument(p, "--metadata", 
                 help = "Select metadata to include in txt file. Multiple labels 
                 can be specified by ',')")
argv = parse_args(p)

all_seurat_integrated = readRDS(file=argv$input)
DefaultAssay(all_seurat_integrated) = argv$assay

metadata = unlist(Map(trimws, strsplit(argv$metadata, ',')))
cat("Saving seurat object as matrix...\n")
seurat_matrix = seurat_obj_to_matrix_txt(all_seurat_integrated, 
                                         assay_type = argv$assay, 
                                         slot_type = argv$slot,
                                         metadata_cols = metadata)
if (is.na(argv$fname)){
  argv$fname = tail(unlist(strsplit(argv$input, '/|.rds')), n=1)
}
write.table(seurat_matrix, file = paste0(argv$output, argv$fname, ".txt"), sep="\t", row.names = TRUE)