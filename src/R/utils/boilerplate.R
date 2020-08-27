# Boilerplate functions for analyzing human / mouse 10X data 
# Author: Koki Sasagawa 

library(Seurat, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(cowplot, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(patchwork)


calculate_QC_metrics = function(seurat_object_list, species){
  if (species == 'human') {
    mito_regex = "^MT-"
    ribo_regex = "^RP[SL]" 
  } else if (species == 'mouse') {
    mito_regex = "^mt-"
    ribo_regex = "^Rp[sl]"
  } else {
    cat("Please specify species = ('human' or 'mouse')")
    return 
  }
  for (i in 1:length(seurat_object_list)) {
    seurat_object_list[[i]][["percent.mt"]] = PercentageFeatureSet(seurat_object_list[[i]], 
                                                                   pattern = mito_regex)
    seurat_object_list[[i]][["percent.ribo"]] = PercentageFeatureSet(seurat_object_list[[i]], 
                                                                     pattern = ribo_regex)
  }
  return(seurat_object_list)
}


calculate_MT_cutoff = function(seurat_object_list, save_directory=NA) {
  percent_mt = data.frame(sample=character(), percent_mt=double())
  for (i in 1:length(seurat_object_list)){
    row_data = data.frame(sample = seurat_object_list[[i]]@meta.data$orig.ident,
                          percent_mt = seurat_object_list[[i]]@meta.data$percent.mt)
    percent_mt = rbind(percent_mt, row_data)
  }
  
  cutoff = get_extreme_outlier(percent_mt$percent_mt)
  
  p1 = ggplot(percent_mt, aes(x=sample, y=percent_mt, fill=sample)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept=cutoff, linetype="dashed", color = "red", size=1)
  print(p1)
  
  p2 = ggplot(percent_mt, aes(x=percent_mt, fill=sample)) +
    geom_histogram() +
    geom_vline(xintercept=cutoff, linetype="dashed", color = "red", size=1)
  print(p2)
  
  if (!is.na(save_directory)){
    pdf(paste0(save_directory, 'percent_mt_boxplot.pdf'), width = 8, height = 5)
    print(p1)
    dev.off()
    
    pdf(paste0(save_directory, 'percent_mt_hist.pdf'), width = 8, height = 5)
    print(p2)
    dev.off()
  }
  
  return(cutoff)
}


get_extreme_outlier = function(values){
  # Extreme outliers are any data values which lie more than 3.0 times the interquartile range 
  upperq = quantile(values)[4]
  lowerq = quantile(values)[2]
  iqr = upperq - lowerq
  cutoff = (iqr * 3) + upperq
  return(cutoff)
}


calculate_cell_cycle_phase = function(seurat_object_list){
  s_genes = cc.genes.updated.2019$s.genes
  g2m_genes = cc.genes.updated.2019$g2m.genes
  for (i in 1:length(seurat_object_list)){
    seurat_object_list[[i]] = CellCycleScoring(seurat_object_list[[i]], 
                                               s.features = s_genes, 
                                               g2m.features = g2m_genes)
  }
  return(seurat_object_list)
}


get_cell_counts = function(seurat_object_list){
  cell_counts = c()
  for (i in 1:length(c(seurat_object_list))){
    cell_counts = c(cell_counts, dim(seurat_object_list[[i]])[2])
  }
  names(cell_counts) = names(seurat_object_list)
  return(cell_counts)
}


percent_reduction = function(original_size, filtered_size, r=2){
  return(round(((original_size - filtered_size) / original_size) * 100, r))
}


add_meta_info = function(seurat_object, meta_info, col_name) {
  meta_list = rep(paste0(meta_info), dim(seurat_object)[2])
  seurat_object = AddMetaData(seurat_object, meta_list, col.name = col_name)
  return(seurat_object)
}


add_multiple_meta_info = function(seurat_object, metadata, cols, fname) {
  # Adds pre-selected metadata information by calling the addMetaInfo function
  # on columns in metadata dataframe
  for (i in cols){
    seurat_object = add_meta_info(seurat_object = seurat_object, 
                                  meta_info = metadata[[i]][metadata$`File Name` == fname],
                                  col_name = sub(' ', '_', i))
  }
  return(seurat_object)
}


get_cell_counts_by_marker = function(seurat_object, marker_list) {
  # Counts the number of cells expressing a particular marker
  counts_vector = vector()
  for (i in marker_list) {
    results = tryCatch(sum(seurat_object@assays$RNA@counts[i,] > 0), 
                       error = function(e) {
                         return(0)})
    counts_vector = append(counts_vector, assign(paste0(i, 'counts'), results))
  }
  names(counts_vector) = marker_list
  return(counts_vector)
}


generate_marker_count_table = function(seurat_object_list,
                                       marker_list, 
                                       save_directory=NA){
  # Generate a table of markers and corresponding number of cells that have some amount of expression > 0 
  for (i in 1:length(seurat_object_list)){
    if (i == 1) {
      counts = get_cell_counts_by_marker(seurat_object_list[[i]], marker_list)
    } else {
      counts = rbind(counts, get_cell_counts_by_marker(seurat_object_list[[i]], marker_list))
    }
  }
  total = get_cell_counts(seurat_object_list) 
  counts = cbind(counts, total)
  rownames(counts) = names(seurat_object_list)
  if (!is.na(save_directory)){
    write.csv(counts,  file = paste0(save_directory, "cell_counts_by_marker.csv")) 
  }
  return(counts)
}


marker_barplots = function(counts, save_directory, y_buffer = 200){
  # Label with '%' symbol
  annotate = function(value) {
    paste(as.character(value)) #, '%') 
  }
  
  for (i in 1:nrow(counts)) {
    # Find percentage of cells expressing the marker 
    percent = round(counts[i,] / counts[i, 'total'], 3) * 100
    percent = lapply(percent, annotate)
    pdf(paste0(save_directory, "cell_count_by_marker_genes_", rownames(counts)[i], ".pdf"))
    par(mar=c(6.5, 4.1, 4.1, 2.1))
    pt2 = barplot(counts[i,], las = 2, 
                  main = paste0("Cell count by marker ", rownames(counts)[i], " (N:", counts[i, 'total'], ')'), 
                  col = rainbow(ncol(counts)), 
                  xpd = FALSE, 
                  ylim = range(pretty(c(0, (counts[i, 'total']+y_buffer)))))
    text(x = pt2, y = counts[i,], label=percent, pos = 3, cex = 0.6, col = "black")
    dev.off()
  }
}


gene_expression_histogram = function(seurat_object,
                                     slot, 
                                     marker_list, 
                                     save_directory,
                                     binsize="Sturges",
                                     rm_zero=F){
  for (i in marker_list){
    # Excluding all 0's, find the quantitles of data 
    values = GetAssayData(seurat_object, slot = slot)[i,]
    if (rm_zero){ 
      values = values[values > 0] 
    }
    quantiles = quantile(values, c(.05, .10, .25, .50, .75, .90, .95))
    pdf(paste0(save_directory, "distribution_of_", i ,"_expression.pdf"))
    h = hist(values, main=paste0(i,' Expression Normalized'),
             xlab=paste0('Normalized Expresson Value Quantiles', '\n',
                         paste(round(quantiles, 3), collapse='  ')),
             breaks=binsize)
    # text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))
    dev.off()
  }
}


subset_by_labels = function(seurat_object, label_value_pairs){
  for (label_val in label_value_pairs){
    label = strsplit(label_val, ':')[[1]][[1]]
    value = strsplit(label_val, ':')[[1]][[2]]
    cat("Subsetting", value, "from", label, "\n")
    seurat_object = subset(seurat_object, subset = !!sym(label) == value)
    cat( "Seurat object dimensions (after subsetting):",
         dim(seurat_object), "\n")
  }
  return(seurat_object)
}


simple_merge = function(seurat_object_list, 
                        normalized = F,
                        n_features = 2000,  
                        regress_out = NULL){
  all_seurat_integrated = merge(x = seurat_object_list[[1]],
                                y = seurat_object_list[2:length(seurat_object_list)],
                                add.cell.ids = as.vector(names(seurat_object_list)),
                                merge.data = normalized,
                                project = "simple_merge")
  if (!normalized){
    all_seurat_integrated = NormalizeData(all_seurat_integrated, 
                                          normalization.method = "LogNormalize", 
                                          scale.factor = 10000, 
                                          verbose = T)
  }
  all_seurat_integrated = FindVariableFeatures(all_seurat_integrated,
                                               selection.method = "vst",
                                               nfeatures = n_features, 
                                               verbose = T)
  all_seurat_integrated = ScaleData(all_seurat_integrated, 
                                    vars.to.regress = regress_out,
                                    verbose = T)
  return(all_seurat_integrated)
}


cca_merge = function(seurat_object_list,
                     n_features = 2000, 
                     k_anchor = 5,
                     k_filter = 200, 
                     k_score = 30, 
                     D = 30, 
                     regress_out = NULL){
                     # references = NULL){
  for (i in 1:length(seurat_object_list)) {
    seurat_object_list[[i]] = NormalizeData(seurat_object_list[[i]], 
                                            normalization.method = "LogNormalize", 
                                            scale.factor = 10000, 
                                            verbose = T)
    seurat_object_list[[i]] = FindVariableFeatures(seurat_object_list[[i]], 
                                                   selection.method = "vst", 
                                                   nfeatures = n_features, 
                                                   verbose = T)
  }
  
  # if (!is.null(references)){
  #   references = which(names(seurat_object_list) %in% references)
  # }
  
  data_anchors = FindIntegrationAnchors(object.list = seurat_object_list,
                                        normalization.method = "LogNormalize",
                                        k.anchor = k_anchor,
                                        k.filter = k_filter,
                                        k.score = k_score,
                                        dims = 1:D,
                                        # reference = references,
                                        verbose = T)
  all_seurat_integrated = IntegrateData(anchorset = data_anchors,
                                        normalization.method = "LogNormalize",
                                        dims = 1:D,
                                        verbose = T)
  DefaultAssay(all_seurat_integrated) = 'integrated'
  all_seurat_integrated = ScaleData(all_seurat_integrated, 
                                    vars.to.regress = regress_out,
                                    verbose = T)
  return(all_seurat_integrated)
}


sctransform_merge = function(seurat_object_list, 
                             n_features = 3000, 
                             keep_var = T, 
                             k_anchor = 5, 
                             k_filter = 200, 
                             k_score = 30, 
                             D = 30, 
                             regress_out = NULL){
                             # references = NULL){
  for (i in 1:length(seurat_object_list)) {
    seurat_object_list[[i]] = SCTransform(seurat_object_list[[i]],
                                          vars.to.regress = regress_out,
                                          variable.features.n = n_features, 
                                          return.only.var.genes = keep_var,
                                          verbose = T)
  }
  # Select Integration features - genes are ranked by the number of datasets they appear in
  data_features = SelectIntegrationFeatures(object.list = seurat_object_list,
                                            nfeatures = n_features,
                                            verbose = T)
  seurat_object_list = PrepSCTIntegration(object.list = seurat_object_list,
                                          anchor.features = data_features, 
                                          verbose = T)
  
  # if (!is.null(references)){
  #   references = which(names(seurat_object_list) %in% references)
  # }
  
  # Identify anchors and integrate datasets
  data_anchors = FindIntegrationAnchors(object.list = seurat_object_list,
                                        normalization.method = "SCT",
                                        anchor.features = data_features,
                                        k.anchor = k_anchor,
                                        k.filter = k_filter,
                                        k.score = k_score,
                                        dims = 1:D,
                                        # reference = references,
                                        verbose = T)
  all_seurat_integrated = IntegrateData(anchorset = data_anchors,
                                        normalization.method = "SCT",
                                        dims = 1:D,
                                        verbose = T)
  DefaultAssay(all_seurat_integrated) = 'integrated'
  return(all_seurat_integrated)
}


########## Figures ##########
standard_scatterplots = function(seurat_object){
  #' Function to generate two scatter plots: count vs percent.mt and count vs nFeature
  plot1 = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plot2 = FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  CombinePlots(plots = list(plot1, plot2), ncol = 2)
}


generate_violin_plots = function(seurat_object_list, 
                                 save_directory=NA, 
                                 show_fig=T){
  plot_list = list()
  for (i in 1:length(seurat_object_list)){
    plot_list[[i]] = VlnPlot(seurat_object_list[[i]], 
                             features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                             ncol = 3, pt.size = 0.5)
    if (show_fig){
      print(plot_list[[i]]) 
    } 
    if (!is.na(save_directory)){
      pdf(paste0(save_directory, levels(seurat_object_list[[i]]), '_violin_plot.pdf'),
          width = 10, height = 7)
      print(plot_list[[i]])
      dev.off()
    }
  }
}


generate_feature_scatter_plots = function(seurat_object_list,
                                          save_directory=NA, 
                                          show_fig=T){
  plot_list = list()
  for (i in 1:length(seurat_object_list)){
    plot_list[[i]] = standard_scatterplots(seurat_object_list[[i]]) 
    if (show_fig){
      print(plot_list[[i]]) 
    } 
    if (!is.na(save_directory)){
      pdf(paste0(save_directory, levels(seurat_object_list[[i]]), '_scatter_plot.pdf'),
          width = 12, height = 7)
      print(plot_list[[i]])
      dev.off()
    }
  }
}

# Following function references this post:
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
modify_vlnplot =  function(seurat_object, 
                           feature, 
                           pt.size = 0, 
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                           ...) {
  p = VlnPlot(seurat_object, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max = function(p){
  ymax = max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot = function(seurat_object, 
                          features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          xaxis_rot = 0,
                          ...) {
  
  plot_list = purrr::map(features,
                         function(x) modify_vlnplot(seurat_object = seurat_object,
                                                    feature = x, ...))
  plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = xaxis_rot), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs = purrr::map_dbl(plot_list, extract_max)
  plot_list = purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p = patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}