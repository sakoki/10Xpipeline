# IO Utility Functions
# Author: Koki Sasagawa 

library(readxl, quietly = TRUE)


create_file_directory = function(file_paths){
  for (i in file_paths) {
    if(!dir.exists(i)) {
      dir.create(i, recursive = TRUE)
    }
  }
}


read_csv_or_xlsx = function(input_file) {
  file_type = tail(unlist(strsplit(input_file, '\\.')), n=1)
  if (file_type == 'csv') {
    return(read.csv(file = input_file))
  } else if (file_type == 'xlsx') {
    return(read_excel(path = input_file))
  } else {
    cat("File must be in the following file formats: 'csv' or 'xlsx'\n")
  }
}


load_10X_files = function(input_directory, save_directory, min_cells = 5, min_features = 200){
  seurat_object_list = c()
  
  cat("Loading data...\n")
  msg = c(paste0("---------- Read data parameters ----------\n", 
                 "Min cell filter: ", min_cells, '\n', 
                 "Min gene filter: ", min_features, '\n'))
  cat(msg)
  
  file_names = c(list.files(input_directory))
  for (file in file_names) {
    seurat_data = Read10X(data.dir = paste0(input_directory, file))
    seurat_object = CreateSeuratObject(counts = seurat_data, 
                                       min.cells = min_cells,
                                       min.features = min_features, 
                                       project = file)
    seurat_object_list = c(seurat_object_list, seurat_object)
    dimensions = seurat_object@assays$RNA@counts@Dim
    msg1 = paste0("Finished loading: ", file, 
                  "\n# of genes: ", dimensions[1], 
                  "\n# of cells: ", dimensions[2], "\n")
    cat(msg1)
    msg = c(msg, msg1)
  }
  
  # Write summary notes 
  cat(msg, file=paste0(save_directory, 'file_summary.txt'), sep='')
  
  names(seurat_object_list) = file_names
  return(seurat_object_list)
} 


load_10X_rds_objects = function(input_directory){
  files = c(list.files(input_directory))
  seurat_object_list = c()
  file_names = c()
  for (i in files){
    tmp = readRDS(paste0(input_directory, i))
    seurat_object_list = c(seurat_object_list, tmp)
    file_names = c(file_names, as.character(tmp@meta.data$orig.ident[1]))
  }
  names(seurat_object_list) = file_names
  return(seurat_object_list)
}

