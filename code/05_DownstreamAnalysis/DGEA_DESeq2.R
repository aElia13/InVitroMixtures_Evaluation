# script with function to run diff exp analysis using DESeq2 R package

# load libraries
library("DESeq2")

# --> following: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

DGEA_DESeq2 <- function(dds,
                        contrast,
                        padj_cutoff = 0.05, # DGEA
                        lfc_cutoff = 0 # DGEA
                        ){

  # run DESeq on dds object
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast)
  
  # differential expression analysis results
  res <- res %>%
    data.frame() %>%
    rownames_to_column("genes") %>%
    dplyr::mutate(contrast = paste(contrast, collapse="_"))
  
  # signifficant diff. expressed genes
  res.filt <- res %>%
    dplyr::filter(padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)
  
  diff_exp <- res.filt %>%
    dplyr::pull("genes")
  
  cutoffs <- list(padj = padj_cutoff,
                  lfc = lfc_cutoff)
  
  # create ranked list of genes from diff exp. by stat column (stat=LFC/se(LFC))
  results_stat <- res %>% 
    dplyr::select("genes", "stat") %>%
    dplyr::arrange(desc(stat))  # arrange by descending order of stat column
  
  rank_list <- results_stat %>% 
    distinct(genes, .keep_all = TRUE) %>% # remove dupl. genes
    na.omit() %>% # remove null values
    deframe() # convert to list where names are gene symbols

  # return results
  return(list(res = res, # results
              res.filt = res.filt, # filtered results
              diff_exp = diff_exp, # differential expressed genes
              rank_list = rank_list, # ranked list of genes for GSEA
              cutoffs = cutoffs, # padj and lfc cutoffs used to filter DEGs
              )
         ) 
}
