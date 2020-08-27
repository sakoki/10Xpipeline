# 10X Analysis Pipeline: Merge Samples
# Author: Koki Sasagawa

library(argparser, quietly = T)
library(Seurat, quietly = T)
library(readxl, quietly = T)
library(here)
source(here::here('src' , 'R', 'utils', 'IO.R'))
source(here::here('src', 'R', 'utils', 'boilerplate.R'))
options(future.globals.maxSize = 100000000000)  # Override memory limit

# Parse command line arguments
p = arg_parser("Merge single cell samples")
p = add_argument(p, "--input", help = "Path to input single cell files")
p = add_argument(p, "--output", help = "Path to save analysis results")
p = add_argument(p, "--tmp", help = "Path to save temporary files")
p = add_argument(p, "--fname", 
                 help = "Name of output seurat object file", 
                 default = "Merged_seurat_obj.rds")
p = add_argument(p, "--subset",
                 help = "Subset the seurat object using an existing label (example: 
                 'cellType:T-cell' will subset by celltype=='T-cell'. Multiple s
                 ubsets can be specified by ',')", 
                 default = 'None')

# Merge parameters
p = add_argument(p, "--merge_type",
                 type = "numeric",
                 help = "Merging logic (1 = 'simple merge', 2 = 'cca', 3 = 'SCTransform')",
                 default = 3)
p = add_argument(p, "--N",
                 type = "numeric", 
                 help = "Number of variable features", 
                 default = 3000)
p = add_argument(p, "--regress_out",
                 help = "Variable to regress out (example: percent.mt. Multiple
                 variables can be specified by ',')", 
                 default = NULL)

# Simple Merge Parameters
p = add_argument(p, "--normalized", 
                 type = "logical", 
                 help = "Merge slot='data' if True, slot='count' if False", 
                 default = F)

# Batch correction parameters
p = add_argument(p, "--keep_var", 
                 type = "logical",
                 help = "Keep only variable genes for SCT normalization step",
                 default = T)
p = add_argument(p, "--D",
                 type = "numeric", 
                 help = "Number of dimensions to use from CCA to specify the 
                 neighbor search space", 
                 default = 30)
# p = add_argument(p, "--reference", 
#                  help = "Sample to use as reference during integration step 
#                  (speed up computation time)", 
#                  default = NULL)
p = add_argument(p, "--k_anchor", 
                 type = "numeric",
                 help = "How many neighbors (k) to use when picking anchors", 
                 default = 5)
p = add_argument(p, "--k_filter", 
                 type = "numeric", 
                 help = "How many neighbors (k) to use when filtering anchors",
                 default = 200)
p = add_argument(p, "--k_score",
                 type = "numeric", 
                 help = "How many neighbors (k) to use when scoring anchors", 
                 default = 30)

# PCA parameters
p = add_argument(p, "--n_pcs",
                 type = "numeric", 
                 help = "Number of principal components to calculate", 
                 default = 50)
argv = parse_args(p)

argv$output = paste0(argv$output, 'merged/')
create_file_directory(argv$output)

# Load Objects
seurat_objects = load_10X_rds_objects(argv$input)

# Set default assay to RNA slot (original data)
for (i in 1:length(seurat_objects)){
  DefaultAssay(seurat_objects[[i]]) = "RNA"
}

if (argv$subset != 'None') {
  argv$subset = unlist(Map(trimws, strsplit(argv$subset, ',')))
  for (i in 1:length(seurat_objects)){
    seurat_objects[[i]] = subset_by_labels(seurat_object = seurat_objects[[i]], 
                                           label_value_pairs = argv$subset)
  }
}

if (!is.null(argv$regress_out)) {
  argv$regress_out = as.vector(mapply(trimws, unlist(strsplit(argv$regress_out, ','))))
}

if (argv$merge_type == 1){  # No Batch Correction (Simple Merge)
  all_seurat_integrated = simple_merge(seurat_objects, 
                                       normalized = argv$normalized, 
                                       n_features = argv$N,
                                       regress_out = argv$regress_out)
  
} else if (argv$merge_type == 2){  # Batch correction (Canonical Correlation Analysis)
  all_seurat_integrated = cca_merge(seurat_objects,
                                    n_features = argv$N,
                                    k_anchor = argv$k_anchor,
                                    k_filter = argv$k_filter,
                                    k_score = argv$k_score, 
                                    D = argv$D,
                                    regress_out = argv$regress_out)
                                    # references = argv$reference)
  
} else if (argv$merge_type == 3){  # Batch correction (Seurat V3 SCTransform)
  all_seurat_integrated = sctransform_merge(seurat_objects,
                                            n_features = argv$N,
                                            keep_var = argv$keep_var,
                                            k_anchor = argv$k_anchor,
                                            k_filter = argv$k_filter,
                                            k_score = argv$k_score, 
                                            D = argv$D,
                                            regress_out = argv$regress_out)
                                            # references = argv$reference)
}

cat("Running PCA...\n")
all_seurat_integrated = RunPCA(all_seurat_integrated, npcs = argv$n_pcs)

# Heuristic to assess how many PC's to use for downstream analysis
pdf(paste0(argv$output, 'elbow_plot.pdf'), width = 10, height = 8)
ElbowPlot(all_seurat_integrated, ndims = argv$n_pcs, reduction = "pca")
dev.off()

# Save the first N PC
write.csv(all_seurat_integrated@reductions$pca@feature.loadings[,1:argv$n_pcs],
          file = paste0(argv$output, paste0('first_', argv$n_pcs, '_PCs.csv')))

cat("Saving Seurat object...\n")
saveRDS(all_seurat_integrated, file=paste0(argv$tmp, argv$fname))

cell_counts = as.data.frame(table(all_seurat_integrated@meta.data$orig.ident))
colnames(cell_counts) = c("orig.ident", "cell_count")
write.csv(cell_counts, file = paste0(argv$output, "cell_count.csv"))