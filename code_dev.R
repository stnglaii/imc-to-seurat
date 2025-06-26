library(tidyverse)
library(Seurat)
library(qs) #archive for saving objects
library(ProfilerAPI2)
options(Seurat.object.assay.version = "v3")

# fetch data
api = ProfilerAPI2::profiler_api(profile = "default")   ## needed for Rmarkdown connection to ProfilerAPI2

# look at data
sample_ids <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgm_20231017_prod_kapusta_0") %>%
  select(sample_id_parent) %>% 
  dplyr::distinct() %>% 
  dplyr::collect()
sample_ids <- sample_ids$sample_id_parent

subject_ids <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgm_20231017_prod_kapusta_0") %>%
  select(subject_id) %>% 
  dplyr::distinct() %>% 
  dplyr::collect()

cohorts <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgm_20231017_prod_kapusta_0") %>%
  select(immucan_cohort) %>% 
  dplyr::distinct() %>% 
  dplyr::collect()

# select sample and ROI
expression <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgm_20231017_prod_kapusta_0") %>%
  filter(sample_id_parent=="IMMU-BC2-1161-FIXT-01") %>% 
  dplyr::distinct() %>% 
  dplyr::collect() %>% 
  mutate(cell_id = paste0("ROI_", image_roi, "_cell_", cell_object_number))

# extract matrix data
matrix <- expression %>%
  select(cell_id, id_imc_analyte, quant_mean_cell_intensity)

# create seurat matrix format
counts <- as.data.frame(pivot_wider(matrix, values_from = quant_mean_cell_intensity, names_from = cell_id))
rownames(counts) <- counts$id_imc_analyte
counts <- counts[,-1]

# create seurat metadata format
metadata <- as.data.frame(expression) %>%
  select(-c(id_immucan_file_name, image_id, id_imc_analyte, quant_mean_cell_intensity, analyte_name, analyte_channel_id, analyte_antibody_target, image_deepcell_channel, analyte_antibody_clone, analyte_antibody_final_concentration)) %>% 
  distinct()
rownames(metadata) <- metadata$cell_id
metadata$cell_id <- NULL

# create seurat object
assay <- CreateAssayObject(counts = counts)
seurat_object <- CreateSeuratObject(counts = assay, assay = "RNA", meta.data = metadata)
seurat_object <- SetIdent(seurat_object, value = "cell_type")

