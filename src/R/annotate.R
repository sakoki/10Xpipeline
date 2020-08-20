# 10X Analysis Pipeline: Annotate cells
# Author: Koki Sasagawa

library(argparser, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(here, quietly = TRUE)
source(here::here('src' , 'R', 'utils', 'IO.R'))
source(here::here('src', 'R', 'utils', 'boilerplate.R'))
options(future.globals.maxSize = 100000000000)   # Override default memory limit on R

# Parse command line arguments
p = arg_parser("Annotate Seurat cells")
p = add_argument(p, "--input", help = "Path to seurat object")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--tmp", help = "Path to save temporary files")
p = add_argument(p, "--annotations", help = "Path to annotation file")
p = add_argument(p, "--fname", 
                 help = "Name of output seurat object file", 
                 default = "Annotated_seurat_obj.rds")
p = add_argument(p, "--subset", 
                 help = "Subset the seurat object using an existing label (example:
                 'cellType:T-cell' will subset by celltype=='T-cell'. Multiple 
                 subsets can be specified by ',')", 
                 default = 'None')
p = add_argument(p, "--annot_name",
                 help = "Name of new annotation (Multiple names can be specified
                 by ',')")
p = add_argument(p, "--column",
                 help = "Column to use for annotation (labels must be in order by
                 cell. Multiple columns can be specified by ',')", 
                 default = 'None')
p = add_argument(p, "--mapping", 
                 help = "Specify mapping of key (col1) to value (col2) for annotation
                 (example: 'col1:col2' where col1 matches a column in the seurat 
                 metadata. Multiple mappings can be specified by ',')", 
                 default = 'None')
p = add_argument(p, "--fig_width",
                 help = "Width of umap figure",
                 type = 'double',
                 default = 12)
p = add_argument(p, "--fig_height",
                 help = "Height of umap figure",
                 type = 'double',
                 default = 7)
argv = parse_args(p)

argv$output = paste0(argv$output, 'annotations/')
create_file_directory(argv$output)

cat("Annotating cells...\n")
annotations = read_csv_or_xlsx(argv$annotations)
all_seurat_integrated = readRDS(file=argv$input)
cat("Seurat object dimensions:", dim(all_seurat_integrated), "\n")

if (argv$subset != 'None') {
  argv$subset = unlist(Map(trimws, strsplit(argv$subset, ',')))
  all_seurat_integrated = subset_by_labels(seurat_object = all_seurat_integrated, 
                                           label_value_pairs = argv$subset)
}

# Annotation
annotation_name = unlist(Map(trimws, strsplit(argv$annot_name, ',')))
if ((argv$mapping != 'None') & (argv$column != 'None')) {
  cat("Invalid: Argument for both mapping and column detected. Can only handle one at a time.\n")
} else if (argv$mapping != 'None') {
  mapping = unlist(Map(trimws, strsplit(argv$mapping, ',')))
  for (i in 1:length(mapping)) {
    col1 = strsplit(mapping[[i]], ':')[[1]][[1]]
    col2 = strsplit(mapping[[i]], ':')[[1]][[2]]
    cat("Unique items in", col1, ":", unique(annotations[[col1]]), "\n")
    cat("Unique items in", col2, ":", unique(annotations[[col2]]), "\n")
    all_seurat_integrated@meta.data[, annotation_name[[i]]] = plyr::mapvalues(x = all_seurat_integrated@meta.data[[col1]],
                                                                              from = annotations[[col1]],
                                                                              to = as.vector(annotations[[col2]]))
  }
} else if (argv$column != 'None') {
  columns = unlist(Map(trimws, strsplit(argv$column, ',')))
  for (i in 1:length(columns)){
    cat("Unique items in", columns[[i]], ":", unique(annotations[[columns[[i]]]]), "\n")
    all_seurat_integrated@meta.data[, annotation_name[[i]]] = as.vector(annotations[[columns[[i]]]])
  }
} else {
  cat("Must specify one annotation type: 'column' or 'mapping'\n")
}

# Figures
for (i in 1:length(annotation_name)){
  pdf(paste0(argv$output, trimws(annotation_name[[i]]), gsub(".rds", "", argv$fname), '_labels.pdf'),
      width = argv$fig_width,
      height = argv$fig_height,
      family = "ArialMT")
  print(DimPlot(all_seurat_integrated, group.by = annotation_name[[i]]))
  dev.off()
}

saveRDS(all_seurat_integrated, file=paste0(argv$tmp, argv$fname))