# script to process output data from deconvolution with DWLS
# and prepare to evaluate by comparing output with ground truth

# Elia Alonso
# 14-10-2025

# load libraries 
library("readr")
library("tibble")
library("dplyr")
library("tidyr")
library("DESeq2")
library("SummarizedExperiment")
library("MultiAssayExperiment")
library("BisqueRNA")

process_output_dwls <- function(res, pval_cutoff){
  # res is a df where rownames are cells and colnames are samples
  # pval_cutoff is the value indicated to build signature matrix
  
  results <- list()  # create empty list to store processed results
  
  est_fractions <- res
  rownames(est_fractions) <- tolower(rownames(est_fractions)) # lowercase for consistency
  print(colnames(est_fractions)) # print samples where fractions were calculated
  
  results[['fractions']] <- SummarizedExperiment(assays = list(cell_fractions = est_fractions))
  
  # arrange in MultiAssayExperiment and return
  results <- MultiAssayExperiment(results,
                                  metadata = list(tool = "DWLS",
                                                  pval.cutoff = pval_cutoff))
  return(results)
}
