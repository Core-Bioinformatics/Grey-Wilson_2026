options(future.globals.maxSize = 8000 * 1024^2)
library(dplyr)
library(ggplot2)
library(Seurat)
library(ClustAssess)
library(harmony)
library(ShinyCell)

ncores <- 14
my_cluster <- parallel::makeCluster(ncores, type = "PSOCK")
doParallel::registerDoParallel(cl = my_cluster)

soMerged <- readRDS("path/to/soMerged.rds")
pc_genes <- readLines("path/to/pc_genes.txt")


counts <- GetAssayData(soMerged, assay = "RNA")
counts <- counts[which(rownames(counts) %in% pc_genes), ]
soMerged <- subset(soMerged, features = rownames(counts))
soMerged <- SCTransform(soMerged, return.only.var.genes = FALSE, verbose = TRUE)

assay_name <- "SCT"
var_features <- soMerged@assays[[assay_name]]@var.features
n_abundant <- 3000
most_abundant_genes <- rownames(soMerged@assays[[assay_name]])[
  order(Matrix::rowSums(soMerged@assays[[assay_name]]), decreasing = TRUE)
][1:n_abundant]


choose_processing_function <- function(proc_type = "default", dt_mtx, npcs = 30, categ = NULL) {
  if (proc_type == "default") {
    return(function(dt_mtx, actual_npcs = 30) {
      actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)
      RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
      embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x
      RhpcBLASctl::blas_set_num_threads(1)
      rownames(embedding) <- rownames(dt_mtx)
      colnames(embedding) <- paste0("PC_", seq_len(ncol(embedding)))
      return(embedding)
    })
  }
  
  if (proc_type == "harmony") {
    return(function(dt_mtx, actual_npcs = 30) {
      actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)
      RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
      embedding <- stats::prcomp(x = dt_mtx, rank. = actual_npcs)$x
      RhpcBLASctl::blas_set_num_threads(1)
      rownames(embedding) <- rownames(dt_mtx)
      colnames(embedding) <- paste0("PC_", seq_len(actual_npcs))
      embedding <- harmony::RunHarmony(embedding, categ, verbose = FALSE)
      return(embedding)
    })
  }
  
  if (proc_type == "cca") {
    return(function(dt_mtx, actual_npcs = 30) {
      actual_npcs <- min(actual_npcs, ncol(dt_mtx) %/% 2)
      dt_mtx <- t(dt_mtx)
      sample_names <- unique(categ)
      
      so_objects <- lapply(sample_names, function(sample_name) {
        ClustAssess::create_seurat_object_default(normalized_expression_matrix = dt_mtx[, categ == sample_name])
      })
      names(so_objects) <- sample_names
      
      unified_so <- merge(x = so_objects[[1]], y = so_objects[-1], add.cell.ids = sample_names)
      unified_so <- NormalizeData(unified_so, verbose = FALSE)
      unified_so <- FindVariableFeatures(unified_so, selection.method = "vst", nfeatures = nrow(dt_mtx), verbose = FALSE)
      unified_so <- ScaleData(unified_so, features = rownames(dt_mtx), verbose = FALSE)
      unified_so@assays$RNA@layers$scale.data <- dt_mtx
      
      RhpcBLASctl::blas_set_num_threads(foreach::getDoParWorkers())
      unified_so <- RunPCA(unified_so, features = rownames(dt_mtx), npcs = actual_npcs, verbose = FALSE, approx = FALSE)
      unified_so <- IntegrateLayers(
        object = unified_so, method = CCAIntegration, orig.reduction = "pca",
        new.reduction = "cca", scale.layer = "scale.data", features = rownames(dt_mtx), verbose = FALSE
      )
      RhpcBLASctl::blas_set_num_threads(1)
      return(unified_so@reductions$cca@cell.embeddings)
    })
  }
}

gene_list <- list(
  "Highly_Variable" = var_features,
  "Most_Abundant" = most_abundant_genes
)

steps_list <- list(
  "Highly_Variable" = c(1500, 1750),
  "Most_Abundant" = c(1500, 1750)
)

expr_matrix <- soMerged@assays$SCT@scale.data
meta <- soMerged@meta.data

RhpcBLASctl::blas_set_num_threads(1)
print('BEGIN ClustAssess')

clustassess_autom <- automatic_stability_assessment(
  expression_matrix = expr_matrix,
  n_repetitions = 30,
  temp_file = "path/to/temp_clustassess.rds",
  n_neigh_sequence = seq(from = 5, to = 50, by = 5),
  resolution_sequence = seq(from = 0.1, to = 2, by = 0.1),
  matrix_processing = choose_processing_function(proc_type = 'harmony', dt_mtx = expr_matrix, categ = meta$Dataset),
  features_sets = gene_list,
  steps = steps_list,
  n_top_configs = 2,
  save_temp = TRUE,
  verbose = TRUE,
  umap_arguments = list(min_dist = 0.3, n_neighbors = 30, metric = "cosine", init = "spectral")
)

saveRDS(clustassess_autom, "path/to/ClustAssess_Output.rds")
parallel::stopCluster(cl = my_cluster)

ClustAssess::write_shiny_app(
  object = soMerged,
  assay_name = "SCT",
  clustassess_object = clustassess_autom,
  project_folder = "path/to/apps/HarmonyCA",
  shiny_app_title = "Merged_Harmony_CA"
)


ca <- readRDS("path/to/ClustAssess_Output.rds")


embedding <- data.frame(ca$Highly_Variable$`1750`$umap)
soMerged@reductions[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(embedding), key = "UMAP_", assay = DefaultAssay(soMerged))

k_values <- c("15", "24", "28", "41")
for (k in k_values) {
  clusters <- ca$Highly_Variable$`1750`$clustering_stability$split_by_k$SLM[[k]]$partitions[[1]]$mb
  ecc <- ca$Highly_Variable$`1750`$clustering_stability$split_by_k$SLM[[k]]$ecc
  
  soMerged@meta.data[[paste0("clusters_", k, "_SLM")]] <- as.factor(clusters)
  soMerged@meta.data[[paste0("ecc_", k)]] <- ecc
}


meta <- soMerged@meta.data
meta$Annotation_Revisions <- 'None'

meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(3,22,12,16)] <- 'Posterior foregut'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(6)] <- 'Proliferative posterior foregut'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(9)] <- 'Proliferative mix'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(11)] <- 'PFG-HE'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(2)] <- 'hepatic endoderm + anterior foregut 2 - need to split'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(10)] <- 'Anterior foregut 2'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(19)] <- 'Anterior foregut 1'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(18)] <- 'Early airway foregut'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(7)] <- 'Anterior foregut 3'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(21)] <- 'Cardiomyocytes'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(14)] <- 'Endothelial'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(5)] <- 'Trachea/Airway progenitor (KRT7,CFTR) + hepatic endoderm'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(1,8,13,4)] <- 'Early HPB'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(15)] <- 'Early HB'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(20,24,17)] <- 'LPM/SHF'
meta$Annotation_Revisions[meta$clusters_24_SLM %in% c(23)] <- 'Mystery cells'

idx <- meta$Annotation_Revisions == 'hepatic endoderm + anterior foregut 2 - need to split'
meta$Annotation_Revisions[idx & meta$clusters_28_SLM %in% c(10,13,15)] <- 'Anterior foregut 2'
meta$Annotation_Revisions[idx & meta$clusters_28_SLM %in% c(1,4,6,12,17)] <- 'Hepatic Endoderm'

idx <- meta$Annotation_Revisions == 'Trachea/Airway progenitor (KRT7,CFTR) + hepatic endoderm'
meta$Annotation_Revisions[idx & meta$clusters_41_SLM %in% c(17)] <- 'Hepatic Endoderm'
meta$Annotation_Revisions[idx & meta$clusters_41_SLM %in% c(19,35)] <- 'Airway foregut'
meta$Annotation_Revisions[idx & meta$clusters_41_SLM %in% c(4,9,10,13,14,20,22,25,33,36)] <- 'Early HPB'

soMerged@meta.data <- meta
saveRDS(soMerged, "path/to/soMerged_Annotated.rds")

scConf <- ShinyCell::createConfig(soMerged)
ShinyCell::makeShinyApp(
  soMerged, scConf, gene.mapping = TRUE, gex.assay = "SCT",
  shiny.dir = "path/to/apps/ShinyCell_Annotated",
  shiny.title = 'Merged_Annotated_ShinyCell'
)
