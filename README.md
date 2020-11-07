# 10X Pipeline

_Author: Koki Sasagawa_

This pipeline is used for loading, processing, and analyzing single cell RNA seq data. Currently, all pipeline components mainly rely on the Seurat R package for handling single cell data.

![pipeline overview](./img/pipeline_overview.png)

# Examples for running scripts:

### Merge.R

```bash
Rscript  ./R/samples_merge.R \
--input "input/" \
--output "output/" \
--tmp "tmp/" \
--fname "name_of_merged_seurat_object.rds" \
--subset "metadata_category:label_to_filter_for" \
--merge_type "3" \
--N "3000" \
--regress_out "metadata_to_regress_out" \
--keep_var "T" \
--D "30" \
--n_pcs "50"
```

### Cluster.R

```bash
Rscript ./R/cluster.R \
--input "./all_sample_cells_merged_seurat_obj.rds" \
--output "./output/" \
--tmp "./tmp/" \
--assay "integrated" \
--N "40" \
--min_dist "0.1" \
--n_epochs "500" \
--lr "1" \
--neighbors "30" \
--res "0.5, 1.0, 1.5" \
--algorithm "4" \
--run_DE "T" \
--assay_DE "SCT" \
--min_FC "0.0" \
--min_percent "0.10" \
--save_seurat "T"
```

### Annotate.R

Example of `annotation.csv`

| seurat_clusters | immune_cell_type |
| :-------------- | :--------------- |
| 0               | T-cell           |
| 1               | B-cell           |
| 2               | Macrophage       |
| ...             | ...              |

```bash
Rscript ./R/annotate.R \
--input "./all_sample_cells_merged_seurat_obj_cluster_res_1_seurat_obj.rds" \
--output "./annotated/" \
--tmp "./tmp/" \
--annotations "./annotations.csv" \
--fname "annotated_all_sample_cells_merged_seurat_obj_cluster_res_1_seurat_obj.rds" \
--annot_name "immune_cell_type" \
--mapping "seurat_clusters:immune_cell_type" \
--fig_width "12" \
--fig_height "7"
```

### Differential_expression.R
```bash
Rscript ./R/differential_expression.R \
--input "./annotated_all_sample_cells_merged_seurat_obj_cluster_res_1_seurat_obj.rds" \
--output "./output/" \
--subset "immune_cell_type:T-cell" \
--assay "SCT" \
--compare_pair "phenotype:Group_A:Group_B,
phenotype:Group_A:Group_C,
phenotype:Group_B:Group_C" \
--min_logFC "0.0" \
--min_percent "0.10" \
--downsample "Inf" \
--only_pos "F" \
--algorithm "wilcox"
```
