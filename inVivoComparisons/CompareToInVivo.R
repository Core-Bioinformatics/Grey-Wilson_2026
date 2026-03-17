options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(pheatmap)

externalDatasets <- list()

output_dir <- "integrated_datasets/Harmony/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

integrate_and_save <- function(dataset_list, output_name) {
  
  if(length(dataset_list) < 2){
    stop("Please provide a list of at least two datasets to integrate.")
  }
  
  soAll <- merge(x = dataset_list[[1]], y = dataset_list[-1])
  
  soAll <- SCTransform(soAll, return.only.var.genes = FALSE, verbose = TRUE)
  DefaultAssay(soAll) <- "SCT"
  soAll <- RunPCA(soAll)
  
  theta <- 0.1
  print(paste0('Working with theta: ', theta))
  
  soAll <- RunHarmony(soAll, group.by.vars = "Dataset", theta = theta, verbose = TRUE)
  soAll <- RunUMAP(soAll, reduction = "harmony", dims = 1:30)
  
  umap_name <- paste0("umap_", format(theta, scientific = FALSE))
  names(soAll@reductions)[names(soAll@reductions) == "umap"] <- umap_name
  
  saveRDS(soAll, file = file.path(output_dir, paste0(output_name, ".rds")))
}
