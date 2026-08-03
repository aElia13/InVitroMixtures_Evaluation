# script to process output data from deconvolution with CIBERSORTx
# and prepare to evaluate by comparing output with ground truth

# Elia Alonso
# 27-08-2025

# load libraries 
library("readr")
library("tibble")
library("dplyr")
library("tidyr")
library("DESeq2")
library("SummarizedExperiment")
library("MultiAssayExperiment")

################################################################################
# CIBERSORTx
################################################################################



process_output_cibersortx <- function(results_folder){
  # results_folder: path to folder that contains output
  
  # init ----
  cx_results <- list() # create empty list to store results
  
  # fractions ----
  cx_fractions <- list.files(path = results_folder,
                             pattern = "_Fractions",
                             full.names = TRUE)
  
  cx_fractions <- read.table(cx_fractions, header = TRUE) %>%
    data.frame() %>%
    tibble::column_to_rownames("Mixture") %>%
    dplyr::select(-c("P.value", "Correlation", "RMSE"))  %>% 
    t()
  
  # cell type names in lowercase for consistency
  rownames(cx_fractions) <- tolower(rownames(cx_fractions))
  
  cx_fractions <- SummarizedExperiment(assays = list(cell_fractions = cx_fractions))
  cx_results[['fractions']] <- cx_fractions
  
  
  # GEPs ----
  cx_geps <- list.files(path = results_folder,
                        pattern = "GEP_GEPs_Filtered.txt",
                        full.names = TRUE)
  
  # load expression matrix
  cx_geps <- read.table(cx_geps, header = TRUE) %>%
    data.frame() %>%
    column_to_rownames("GeneSymbol")
  
  # cell type names in lowercase for consistency
  colnames(cx_fractions) <- tolower(colnames(cx_fractions))
  
  # arrange in SE object and return
  cx_results[['geps']] <- SummarizedExperiment(assays = list(counts = cx_geps))
  cx_results <- MultiAssayExperiment(cx_results,
                                     metadata = list(tool = "CIBERSORTx"))
  
  return(cx_results)
}
