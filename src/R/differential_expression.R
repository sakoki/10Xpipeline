# 10X Analysis Pipeline: Differential Expression
# Author: Koki Sasagawa

library(argparser, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(limma, quietly = TRUE)
library(here, quietly = TRUE)
source(here::here('src' , 'R', 'utils', 'IO.R'))
source(here::here('src', 'R', 'utils', 'boilerplate.R'))
options(future.globals.maxSize= 100000000000)  # Override memory limit

filter_genes = function(de_results){
  if (length(grep('^RP[SL]', rownames(de_results))) > 0){  # Remove ribosomal genes
    de_results = de_results[-grep('^RP[SL]', rownames(de_results)), ]
  }
  # Ashurst et al., 2005: “If an approved symbol is not available for a gene locus, 
  # an interim internal identifier is used, which is usually in the format 
  # clonename.number, e.g. RP11-694B14.5”.
  if (length(grep('^RP11-', rownames(de_results))) > 0){  # Remove unnamed genes
    de_results = de_results[-grep('^RP11-', rownames(de_results)), ]
  }
  if (length(grep('^MT-', rownames(de_results))) > 0){  # Remove mitochondrial processing factors
    de_results = de_results[-grep('^MT-', rownames(de_results)), ]
  }
  return(de_results)
}

# Parse command line arguments
p = arg_parser("Differential Expression 10X samples")
p = add_argument(p, "--input", help = "Path to Seurat object")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--subset",
                 help = "Subset the seurat object using an existing label (example: 
                 'cellType:T-cell' will subset by celltype=='T-cell'. Multiple s
                 ubsets can be specified by ',')", 
                 default = 'None')
p = add_argument(p, "--assay", help = "Assay to use data for DE", default = "RNA")
p = add_argument(p, "--compare_pair", 
                 help = "Within a category, specify group1 to compare against group2 
                 (example: 'category:group1:group2'. Multiple pairings can be specified
                 by ',')", 
                 default = "None")
p = add_argument(p, "--compare_all",
                 help = "Perform all pairwise comparisons of groups within a specified
                 category (Multiple categories can be specified by ',')", 
                 default = "None")

# DE parameters
p = add_argument(p, "--min_logFC", 
                 help = "Minimum logFC for a differentially expressed gene", 
                 type = 'double',
                 default = 0.25)
p = add_argument(p, "--min_percent",
                 help = "Minimum percentage of cells to be expressing a particular gene",
                 type = 'double', 
                 default = 0.10)
p = add_argument(p, "--downsample",
                 help = "Downsample each identity class", 
                 type = 'double',
                 default = Inf)
p = add_argument(p, "--only_pos", 
                 help = "If TRUE, keep only positive or upregulated genes", 
                 type = 'logical')
p = add_argument(p, "--algorithm",
                 help = "Algorithm for DE analysis (currently supporting those
                 supported by Seurat)")
argv = parse_args(p)
argv$output = paste0(argv$output, 'DE/')
create_file_directory(argv$output)

cat("Running DE...\n")
all_seurat_integrated = readRDS(file=argv$input)
cat("Seurat object dimensions:", dim(all_seurat_integrated), "\n")

if (argv$subset != 'None') {
  argv$subset = unlist(Map(trimws, strsplit(argv$subset, ',')))
  all_seurat_integrated = subset_by_labels(seurat_object = all_seurat_integrated, 
                                           label_value_pairs = argv$subset)
}

DefaultAssay(all_seurat_integrated) = argv$assay
cat(paste0("Selected assay: ", DefaultAssay(all_seurat_integrated), "\n"))
cat(paste0("Selected algorithm: ", argv$algorithm, "\n"))

# DE_analysis
if (argv$compare_pair != "None"){
  comparisons = unlist(Map(trimws, strsplit(argv$compare_pair, ',')))
  for (i in 1:length(comparisons)) {
    category = strsplit(comparisons[[i]], ':')[[1]][[1]]
    group1 = strsplit(comparisons[[i]], ':')[[1]][[2]]
    group2 = strsplit(comparisons[[i]], ':')[[1]][[3]]
    cat(paste0('category: ', category, ' --> compare ', group1, ' against ', group2, '\n'))
    all_seurat_integrated = SetIdent(all_seurat_integrated, value = category)
    all_pheno_markers = FindMarkers(all_seurat_integrated,
                                    ident.1 = group1,
                                    ident.2 = group2,
                                    logfc.threshold = argv$min_logFC, 
                                    min.pct = argv$min_percent, 
                                    max.cells.per.ident = argv$downsample, 
                                    only.pos = argv$only_pos,
                                    test.use = argv$algorithm,
                                    verbose=T)
    all_pheno_markers = as.data.frame(all_pheno_markers)
  
  # Median center and sort in descending order
  all_pheno_markers$avg_logFC_median_centered = all_pheno_markers$avg_logFC - median(all_pheno_markers$avg_logFC)
  all_pheno_markers = all_pheno_markers %>% arrange(-avg_logFC_median_centered)
  write.csv(all_pheno_markers, file = paste0(argv$output, 
                                             group1, '_against_', group2, 
                                             '_', category, 
                                             '_', argv$algorithm,
                                             '_diff_exp_genes.csv'))
  all_pheno_markers = filter_genes(de_results = all_pheno_markers)
  write.csv(all_pheno_markers, file = paste0(argv$output, 
                                             group1, '_against_', group2,
                                             '_', category, 
                                             '_', argv$algorithm, 
                                             '_diff_exp_genes_filtered.csv'))
  pdf(paste0(argv$output, 
             group1, '_against_', group2,
             '_', category, 
             '_', argv$algorithm, 
             "_top_genes.pdf"),
      width = 3,
      height = 15,
      family = "ArialMT")
  print(StackedVlnPlot(seurat_object = subset(all_seurat_integrated, 
                                              subset = (!!sym(category) == group1) | 
                                                (!!sym(category) == group2)),
                       features = rownames(head(all_pheno_markers, 15)),
                       xaxis_rot = 45))
  dev.off()
  }
}

if (argv$compare_all != "None"){
  comparisons = unlist(Map(trimws, strsplit(argv$compare_all, ',')))
  for (i in 1:length(comparisons)) {
    category = comparisons[[i]]
    cat(paste0('category: ', category, '\n'))
    all_seurat_integrated = SetIdent(all_seurat_integrated, value = category)
    all_pheno_markers = FindAllMarkers(all_seurat_integrated,
                                       logfc.threshold = argv$min_logFC, 
                                       min.pct = argv$min_percent, 
                                       max.cells.per.ident = argv$downsample, 
                                       only.pos = argv$only_pos,
                                       test.use = argv$algorithm,
                                       verbose=T)
    all_pheno_markers = as.data.frame(all_pheno_markers)
    write.csv(all_pheno_markers, file = paste0(argv$output, 
                                               category, 
                                               '_', argv$algorithm,
                                               '_diff_exp_genes.csv'))
    all_pheno_markers = filter_genes(de_results = all_pheno_markers)
    write.csv(all_pheno_markers, file = paste0(argv$output, 
                                               category, 
                                               '_', argv$algorithm, 
                                               '_diff_exp_genes_filtered.csv'))
    pdf(paste0(argv$output, 
               category, 
               '_', argv$algorithm, 
               "_top_genes.pdf"),
        width = dim(unique(all_seurat_integrated[[category]]))[1],  # dynamically update width
        height = 15,
        family = "ArialMT")
    print(StackedVlnPlot(seurat_object = all_seurat_integrated, 
                         features = head(all_pheno_markers, 15)$gene, 
                         xaxis_rot = 45))
    dev.off()
  }
}

