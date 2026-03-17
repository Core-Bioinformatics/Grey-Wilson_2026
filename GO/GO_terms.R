library(reshape2)
library(viridis)
library(Seurat)
library(gprofiler2)
library(ggplot2)
library(reshape2)
library(viridis)
library(RColorBrewer)
so <- readRDS('object.rds')
head(so@meta.data)
Idents(so) <- so@meta.data$Annotation_Revisions
table(so@meta.data$Annotation_Revisions)
Idents(so) <- so@meta.data$Annotation_Revisions
lfc <- 2.5

so <- PrepSCTFindMarkers(so)
CMMarkers <- FindMarkers(so,ident.1 = 'Cardiomyocytes', logfc.threshold=lfc)
sigCM <- rownames(CMMarkers)[CMMarkers$p_val_adj<0.05]
EndMarkers <- FindMarkers(so,ident.1 = 'Endothelial', logfc.threshold=lfc)
sigEnd <- rownames(EndMarkers)[EndMarkers$p_val_adj<0.05]
allGenes <- rownames(so)
resultsCM <- gost(query = sigCM,
                  organism = "hsapiens",
                  significant = TRUE,
                  user_threshold = 0.05,
                  correction_method = "fdr",
                  custom_bg = allGenes)$result
resultsCM <- apply(resultsCM,2,as.character)
write.table(resultsCM, file = '2D-CM_GEO.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
resultsEnd <- gost(query = sigEnd,
                   organism = "hsapiens",
                   significant = TRUE,
                   user_threshold = 0.05,
                   correction_method = "fdr",
                   custom_bg = allGenes)$result
resultsEnd <- apply(resultsEnd,2,as.character)
write.table(resultsEnd, file = '2D-End_GEO.tsv', quote = FALSE, row.names = FALSE, sep = '\t')

lsMarkers <- FindMarkers(so,ident.1 = 'LPM/SHF', logfc.threshold=lfc)
sigLS <- rownames(lsMarkers)[lsMarkers$p_val_adj<0.05]
resultsLS <- gost(query = sigLS,
                   organism = "hsapiens",
                   significant = TRUE,
                   user_threshold = 0.05,
                   correction_method = "fdr",
                   custom_bg = allGenes)$result
resultsLS <- apply(resultsLS,2,as.character)

eHPBMarkers <- FindMarkers(so,ident.1 = 'Early HPB', logfc.threshold=lfc)
sig_eHPB <- rownames(eHPBMarkers)[eHPBMarkers$p_val_adj<0.05]
resultseHPB <- gost(query = sig_eHPB,
                  organism = "hsapiens",
                  significant = TRUE,
                  user_threshold = 0.05,
                  correction_method = "fdr",
                  custom_bg = allGenes)$result
resultseHPB <- apply(resultseHPB,2,as.character)


HEMarkers <- FindMarkers(so,ident.1 = 'Hepatic Endoderm', logfc.threshold=lfc)
sigHE <- rownames(HEMarkers)[HEMarkers$p_val_adj<0.05]
resultsHE <- gost(query = sigHE,
                  organism = "hsapiens",
                  significant = TRUE,
                  user_threshold = 0.05,
                  correction_method = "fdr",
                  custom_bg = allGenes)$result
resultsHE <- apply(resultsHE,2,as.character)

resultsLS <- as.data.frame(resultsLS)
resultsCM <- as.data.frame(resultsCM)
resultsEnd <- as.data.frame(resultsEnd)
resultseHPB <- as.data.frame(resultseHPB)
resultsHE <- as.data.frame(resultsHE)

terms <- c('GO:0072359','GO:0001944','GO:0035295','GO:0001568','GO:0001525','GO:0043542','GO:0002040','GO:0045446',
           'GO:0007507','GO:0060047','GO:0048738','GO:0014706','GO:0045214','GO:0030016','GO:0030017','GO:0043292','GO:0016459','WP:WP383')
resultsCMfilter <- resultsCM[resultsCM$term_id %in% terms,]
resultsEndfilter <- resultsEnd[resultsEnd$term_id %in% terms,]
resultsLSfilter <- resultsLS[resultsLS$term_id %in% terms,]
resultseHPBfilter <- resultseHPB[resultseHPB$term_id %in% terms,]
resultsHEfilter <- resultsHE[resultsHE$term_id %in% terms,]

resultsCMfilter$CellType <- 'Cardiomyocyte'
resultsEndfilter$CellType <- 'Endothelial'
resultsLSfilter$CellType <- 'LPM/SHF'
resultseHPBfilter$CellType <- 'Early HPB'
resultsHEfilter$CellType <- 'Hepatic Endoderm'

enrichment_all <- rbind(resultsCMfilter, resultsEndfilter,resultsLSfilter,resultseHPBfilter,resultsHEfilter)
enrichment_all$term_name <- factor(enrichment_all$term_name, levels = unique(enrichment_all$term_name))
enrichment_all$enrichment_ratio <- as.numeric(enrichment_all$intersection_size) / as.numeric(enrichment_all$term_size)
enrichment_all$log_pval <- -log10(as.numeric(enrichment_all$p_value))
enrichment_all$CellType <- factor(enrichment_all$CellType,
                                  levels = c("Endothelial", "Cardiomyocyte",'LPM/SHF','Early HPB','Hepatic Endoderm'))
enrichment_all$term_name <- factor(enrichment_all$term_name,
                                   levels = unique(enrichment_all$term_name))

enrichment_all <- enrichment_all[enrichment_all$enrichment_ratio>0.1,]

ggplot(enrichment_all, aes(x = CellType,
                           y = term_name,
                           size = enrichment_ratio,
                           color = log_pval)) +
  geom_point(alpha = 0.8) +
  scale_color_distiller(palette = "YlGnBu", direction = 1, name = "-log10(p-value)") +
  scale_size(range = c(3, 10), name = "Enrichment ratio") +
  labs(
    x = "Cell Type",
    y = "GO Term",
    title = "GO Enrichment Comparison: Endothelial vs Cardiomyocytes"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

library(reshape2)
library(viridis)
library(Seurat)
library(gprofiler2)
library(ggplot2)
library(RColorBrewer)

so <- readRDS('object.rds')
Idents(so) <- so@meta.data$Annotation_Revisions
lfc <- 2.5

so <- PrepSCTFindMarkers(so)
allGenes <- rownames(so)

cell_types <- unique(so@meta.data$Annotation_Revisions)
terms <- c('GO:0072359','GO:0001944','GO:0035295','GO:0001568','GO:0001525','GO:0043542',
           'GO:0002040','GO:0045446','GO:0007507','GO:0060047','GO:0048738','GO:0014706',
           'GO:0045214','GO:0030016','GO:0030017','GO:0043292','GO:0016459','WP:WP383')

enrichment_list <- list()
for(ct in cell_types) {
  markers <- FindMarkers(so, ident.1 = ct, logfc.threshold = lfc)
  sig_genes <- rownames(markers)[markers$p_val_adj < 0.05]
  
  results <- gost(query = sig_genes, organism = "hsapiens", significant = TRUE,
                  user_threshold = 0.05, correction_method = "fdr", custom_bg = allGenes)$result
  
  if(!is.null(results)) {
    results <- as.data.frame(apply(results, 2, as.character))
    results_filter <- results[results$term_id %in% terms, ]
    if(nrow(results_filter) > 0) {
      results_filter$CellType <- ct
      enrichment_list[[ct]] <- results_filter
    }
  }
}

enrichment_all <- do.call(rbind, enrichment_list)

all_terms <- unique(enrichment_all$term_name)
complete_grid <- expand.grid(CellType = cell_types, term_name = all_terms, stringsAsFactors = FALSE)
enrichment_all <- merge(complete_grid, enrichment_all, by = c("CellType", "term_name"), all.x = TRUE)

enrichment_all$enrichment_ratio <- as.numeric(enrichment_all$intersection_size) / as.numeric(enrichment_all$term_size)
enrichment_all$log_pval <- -log10(as.numeric(enrichment_all$p_value))

enrichment_all$CellType <- factor(enrichment_all$CellType, levels = cell_types)
enrichment_all$term_name <- factor(enrichment_all$term_name, levels = all_terms)


term_order <- c("myofibril", "contractile muscle fiber",'Striated muscle contraction pathway','sarcomere organization',
                       'striated muscle tissue development','cardiac muscle tissue development','heart development',"myosin complex",
                       'sarcomere','heart contraction','circulatory system development','endothelial cell differentiation','endothelial cell migration',
                       'sprouting angiogenesis','tube development','vasculature development','angiogenesis','blood vessel development'
                       )
enrichment_all$term_name <- factor(enrichment_all$term_name, levels = term_order)

p <- ggplot(enrichment_all, aes(x = CellType, y = term_name, size = enrichment_ratio, color = log_pval)) +
  geom_point(alpha = 0.8) +
  scale_color_distiller(palette = "YlGnBu", direction = 1, name = "-log10(p-value)", na.value = "grey90") +
  scale_size(range = c(3, 10), name = "Enrichment ratio") +
  labs(x = "Cell Type", y = "GO Term", title = "GO Enrichment Comparison") +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())
ggsave("Plots/LGW_Revisions_GOterms.svg", plot = p, device = "svg", width = 10, height = 7)
