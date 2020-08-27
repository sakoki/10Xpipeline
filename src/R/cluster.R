# 10X Analysis Pipeline: Cluster cells 
# Author: Koki Sasagawa 

library(argparser, quietly = TRUE)
library(Seurat, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(here, quietly = TRUE)
source(here::here('src' , 'R', 'utils', 'IO.R'))
options(future.globals.maxSize= 100000000000)  # Override memory limit

# Parse command line arguments
p = arg_parser("Cluster cells")
p = add_argument(p, "--input", help = "Path to seurat object")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--tmp", help = "Path to save temporary files")

# Parameters for UMAP
p = add_argument(p, "--assay", help = "Default assay")
p = add_argument(p, "--N", 
                 type = "numeric", 
                 help = "N dimensions to use for UMAP",
                 default = 30)
p = add_argument(p, "--min_dist", 
                 type = "numeric", 
                 help = "Parameter to adjustclustering size",
                 default = 0.3)
p = add_argument(p, "--n_epochs", 
                 type = "numeric",
                 help = "Number of epochs (UMAP)", 
                 default = 500)
p = add_argument(p, "--lr",
                 type = "numeric",
                 help = "Learning rate (UMAP)", 
                 default = 1)
p = add_argument(p, "--neighbors",
                 type = "numeric", 
                 help = "Number of neighboring points used in local approximations
                 of manifold structure", 
                 default = 30)

# Parameters for graph clustering 
p = add_argument(p, "--res",
                 help = "Resolution parameter for graph partitioning algorithm 
                 (Multiple values can be specified by ',')")
p = add_argument(p, "--algorithm",
                 type = "numeric",
                 help = "Graph partitioning algorithm",
                 default = 4)

# Parameters for DE
p = add_argument(p, "--run_DE",
                 type = 'logical',
                 help = "Run Seurats default differential expression between meta clusters",
                 default = F)
p = add_argument(p, "--assay_DE",
                 help = "Assay to use for DE",
                 default = 'SCT')
p = add_argument(p, "--min_FC",
                 type = 'numeric',
                 help = "Minimum fold change between clusters",
                 default = 0.10)
p = add_argument(p, "--min_percent",
                type = 'numeric',
                help = "Minimum percent of expression",
                default = 0.10)
# Save files
p = add_argument(p, "--save_seurat",
                 type = 'logical',
                 help = "Save seurat object for later",
                 default = F)
argv = parse_args(p)

# Configure directories
argv$output = paste0(argv$output, 'clustered/')
create_file_directory(argv$output)

# Load file
all_seurat_integrated = readRDS(file=argv$input)
DefaultAssay(all_seurat_integrated) = argv$assay

# UMAP
cat("Running UMAP...\n")
all_seurat_integrated = RunUMAP(all_seurat_integrated,
                                min.dist = argv$min_dist,
                                n.epochs = argv$n_epochs,
                                learning.rate = argv$lr,
                                n.neighbors = argv$neighbors,
                                dims=1:argv$N)
argv$min_dist = gsub('\\.', '-', as.name(argv$min_dist))
pdf(paste0(argv$output,
           'umap_mdist', argv$min_dist,
           '_neighbors', argv$neighbors,
           '_epoch', argv$n_epochs,
           '_lr', argv$lr,
           '_dims', argv$N,
           '.pdf'), width = 14, height = 7)
DimPlot(all_seurat_integrated, reduction = "umap", group.by = "orig.ident")
dev.off()

# Graph Clustering
cat("Building graph...\n")
all_seurat_integrated = FindNeighbors(all_seurat_integrated,
                                      dims = 1:argv$N,
                                      verbose = TRUE)

# Adjust resolution by number of cells --> larger resolution detects more subclusters 
res_list = unlist(Map(trimws, strsplit(argv$res, ',')))
for (i in 1:length(res_list)) {
  res = as.numeric(res_list[i])
  subfolder = paste0(argv$output, "res", res, "/")
  create_file_directory(subfolder)
  
  cat(paste0('Clustering with leiden algorithm at resolution ', res, '\n'))
  clustered = FindClusters(all_seurat_integrated, 
                           resolution = res, 
                           algorithm = argv$algorithm)
  res = gsub('\\.', '-', as.name(res))
  
  # Figures
  pdf(paste0(subfolder, 'clustered_umap_leiden', res, '.pdf'), width = 11, height = 7)
  print(DimPlot(clustered, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5))
  dev.off()

  pdf(paste0(subfolder, 'ind_seurat_clusters_umap_res_', res, '.pdf'), width=15, height=12)
  print(DimPlot(clustered, split.by = "seurat_clusters", ncol=5))
  dev.off()
  
  # DE Markers of clusters
  if (argv$run_DE) {
    DefaultAssay(clustered) = argv$assay_DE
    cat("Finding DE markers between clusters...\n")
    cell_markers = FindAllMarkers(clustered, 
                                  only.pos = TRUE, 
                                  min.pct = argv$min_percent, 
                                  logfc.threshold = argv$min_FC)
    write.csv(cell_markers, file = paste0(subfolder, 'cluster_markers_leiden', res, '.csv'))
    N = 5
    top_makers_by_cluster = cell_markers %>% group_by(cluster) %>% top_n(n=N, wt = avg_logFC)
    write.csv(top_makers_by_cluster, file = paste0(subfolder, "clusters_top_", N, "_markers_leiden", res, ".csv"))
  }
  
  if (argv$save_seurat) {
    fname = tail(unlist(strsplit(argv$input, '/|.rds')), n=1)
    saveRDS(clustered, file=paste0(argv$tmp, fname, "_cluster_res_", res, "_seurat_obj.rds"))
  }
}
