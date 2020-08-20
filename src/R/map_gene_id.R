# Mapper for converting gene nomenclature (human or mouse)
# Author: Koki Sasagawa
library(argparser, quietly = TRUE)
library(here)
source(here::here('src' , 'R', 'utils', 'IO.R'))

# Parse command line arguments
p = arg_parser("Differential Expression 10X samples")
p = add_argument(p, "--input", help = "Path to files")
p = add_argument(p, "--output", help = "Path to save results")
p = add_argument(p, "--col_index", 
                 type = "double",
                 help = "Column index of gene names (Indexing start from 1)",
                 default = 1)
p = add_argument(p, "--species",
                 help = "Specify the species for gene names")
p = add_argument(p, "--mapping_key",
                 help = "Mapping gene name from ENSEMBL, ENTREZID, SYMBOL")
p = add_argument(p, "--mapping_cols", 
                 help = "Mapping gene names to ENSEMBL, ENTREZID, SYMBOL (Multiple
                 formats can be specified by ',')")
argv = parse_args(p)

map_ID = function(df, reference, mapping_key, mapping_col, col=1, save_directory=NA, 
                  fname=NA, verbose=T){
  df[[mapping_key]] = as.character(df[[col]])
  
  if (mapping_key == "ENSEMBL"){  # Remove version number which is sometimes present
    df[[mapping_key]] = gsub("\\.[0-9]*$", "", df[[mapping_key]])
  }

  msg = c(paste0('Number of gene symbols: ', length(df[[mapping_key]]), '\n',
                 'Number of unique gene symbols: ', length(unique(df[[mapping_key]])), '\n'))
  
  for (cols in mapping_col) {
    df[[cols]] = mapIds(reference, 
                        df[[mapping_key]],
                        keytype = mapping_key,
                        column = cols,
                        multiVals = 'first')
    
    msg = c(msg, paste0(cols, ' unsuccessful mappings: ', sum(is.na(df[[cols]])), '\n'))
  }
  
  if (verbose) {
    cat(msg)
  }
  
  if (!is.na(save_directory)) {
    cat(msg, file=paste0(save_directory, fname, '_mapping_results.txt'), sep='') 
  }
  
  return(df)
}

# Load Data
df = read_csv_or_xlsx(input_file = argv$input)

fname = tail(unlist(strsplit(argv$input, '/|.csv')), n=1)
mapping_col = unlist(Map(trimws, strsplit(argv$mapping_col, ',')))

if (grepl("[Hh]uman", argv$species)){
  library("org.Hs.eg.db", quietly = T)
  mapping = map_ID(df, 
                   reference = org.Hs.eg.db,
                   mapping_key = argv$mapping_key,
                   mapping_col = mapping_col,
                   col=1, 
                   save_directory = argv$output,
                   fname = fname,
                   verbose = T)
  write.csv(mapping, file = paste0(argv$output, fname, "mapped.csv"))
}

if (grepl("[Mm]ouse", argv$species)){
  library("org.Mm.eg.db", quietly = T)
  mapping = map_ID(df, 
                   reference = org.Mm.eg.db,
                   mapping_key = argv$mapping_key,
                   mapping_col = mapping_col,
                   col=1, 
                   save_directory = argv$output,
                   fname = fname,
                   verbose = T)
  write.csv(mapping, file = paste0(argv$output, fname, "mapped.csv"))
}
