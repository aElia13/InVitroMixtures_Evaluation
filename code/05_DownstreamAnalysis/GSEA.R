# Script to impolement functions to run gene set enrichment analysis
# on pre-renaked list of genes based on dfferential gene expression analysis

# Elia Alonso
# 22.01.2026

# Load libraries
library("clusterProfiler")
library("fgsea")

# Define function:
GSEA <- function(rank_list,
                 t2g, 
                 padj_cutoff = 0.1,
                 title = NULL){
  
  results <- list() # init empty list to store results
  
  # run GSEA
  gsea_res <- clusterProfiler::GSEA(geneList = rank_list,
                                    TERM2GENE = t2g,
                                    by = "fgsea",
                                    eps = 0, # You can set the `eps` argument to zero for better estimation.
                                    pvalueCutoff = padj_cutoff,
                                    # minGSSize = min_gs_size,
                                    # maxGSSize = max_gs_size,
                                    pAdjustMethod = "BH") 
  
  gsea_res <- gsea_res@result %>% 
    data.frame() %>% 
    dplyr::mutate(source = title)
  
  results[["gsea_res"]] <- gsea_res
  results[["gsea_cutoff"]] <- list(padj = padj_cutoff)
  
  return(results)
}
