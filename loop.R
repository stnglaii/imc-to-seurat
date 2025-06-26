library(tidyverse)
library(Seurat)
library(qs) #archive for saving objects
library(ProfilerAPI2)
options(Seurat.object.assay.version = "v3")

# fetch data
api = ProfilerAPI2::profiler_api(profile = "default")   ## needed for Rmarkdown connection to ProfilerAPI2

# glimpse at data
sample_ids <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgmpanel1_20250313_prod_pruss_1") %>%
  select(sample_id_parent) %>% 
  dplyr::distinct() %>% 
  dplyr::collect()

# assess which samples have been converted already
sample_ids <- sample_ids$sample_id_parent
if (!dir.exists("seurat_files")) {dir.create("seurat_files")}
done <- c(done, tools::file_path_sans_ext(basename(list.files(path = "~/ImmuCan/seurat_files"))))
qsave(done, "./finished_samples.qs")
sample_ids <- setdiff(sample_ids, done)
qsave(sample_ids, "./samples_left.qs")

# initate empty vector to save samples that throw errors
sample_errors <- c()

# conversion loop
for (sample in sample_ids) {
  result <- try({
    # select sample and ROI
    expression <- dplyr::tbl(api$conn,"imi2_immucan_img_imc_csgmpanel1_20250313_prod_pruss_1") %>%
      filter(sample_id_parent==sample) %>% 
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
    
    ROIs <- sort(unique(seurat_object$image_roi))
    
    # create FOV with coordinates for each ROI
    for (current_roi in ROIs) {
      coords <- FetchData(seurat_object, vars = c("image_roi", "cell_pos_x", "cell_pos_y")) %>% 
        filter(image_roi == current_roi) %>% 
        select(-image_roi)
      seurat_object@images[[paste0("ROI_", current_roi)]] <- CreateFOV(coords, type = "centroids", assay = "RNA")
    }
    
    # save seurat object and upload to datalake
    qsave(seurat_object, paste0("~/ImmuCan/seurat_files/",sample, ".qs"))
    upload_file_id <- api$data_lake$upload_file("folder-49d71b87-0dfe-4a4c-bf86-40d25aa8ed46",paste0("~/ImmuCan/seurat_files/", sample, ".qs"))
    print(paste("Done with", sample))
  })
  
  # save samples that throw errors
  if (inherits(result, "try-error")) {
    sample_errors <- c(sample_errors, sample)
    next  # skip to the next iteration if an error occurred
  }
}
