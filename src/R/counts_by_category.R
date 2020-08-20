# 10X Analysis Pipeline: Get cell counts by category
# Author: Koki Sasagawa 

library(argparser, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
source(here::here('src' , 'R', 'utils', 'IO.R'))
source(here::here('src' , 'R', 'utils', 'boilerplate.R'))
options(future.globals.maxSize= 100000000000)   # Override default memory limit on R

# Parse commandline arguments
p = arg_parser("Split the data into categories by one label, then count the number
               of cells belonging to categories of another label")
p = add_argument(p, "--input", help = "Path to seurat object")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--subset", 
                 help = "Subset the seurat object using an existing label (example:
                 'cellType:T-cell' will subset by celltype=='T-cell'. Multiple 
                 subsets can be specified by ',')", 
                 default = 'None')
p = add_argument(p, "--label_pairs", 
                 help = "Specify label to split object by, then label for counting
                 (example: 'label1:label2'. Multiple label pairs can be specified 
                 by ',')")
argv = parse_args(p)

argv$output = paste0(argv$output, 'cell_counts/')
create_file_directory(argv$output)

all_seurat_integrated = readRDS(file=argv$input)
cat("Seurat object dimensions:", dim(all_seurat_integrated), "\n")

if (argv$subset != 'None') {
  argv$subset = unlist(Map(trimws, strsplit(argv$subset, ',')))
  all_seurat_integrated = subset_by_labels(seurat_object = all_seurat_integrated, 
                                           label_value_pairs = argv$subset)
}

label_pairs = unlist(Map(trimws, strsplit(argv$label_pairs, ',')))
for (i in 1:length(label_pairs)) {
  label1 = strsplit(label_pairs[[i]], ':')[[1]][[1]]
  label2 = strsplit(label_pairs[[i]], ':')[[1]][[2]]
  cat("Getting cell counts by label pairs:", label1, "&", label2, "\n")
  counts = c()
  sub_names = c()
  for (i in SplitObject(all_seurat_integrated, split.by = label1)) {
    sub_names = c(sub_names, unique(i@meta.data[, label1]))
    counts = rbind(counts, table(i[[label2]]))
  }
  rownames(counts) = sub_names
  write.csv(counts, file = paste0(argv$output, label1, '_', label2, "_cell_type_counts.csv"))
}
