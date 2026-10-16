# script to process output data from deconvolution with BayesPrism
# and prepare to evaluate by comparing output with ground truth

# Elia Alonso
# 25-07-2025

# load libraries 
library("readr")
library("tibble")
library("dplyr")
library("tidyr")
library("BayesPrism")
library("InstaPrism")
library("DESeq2")
library("SummarizedExperiment")
library("MultiAssayExperiment")

process_output_bayesprism <- function(bp.res){
  key_cell <- bp.res@prism@key
  bp_results <- list()  # create empty list to store processed results
  
  # fractions ----
  # posterior mean of cell type fraction theta
  bp_fractions <- get.fraction(bp = bp.res,
                               which.theta = "final",
                               state.or.type = "type") %>%
    t() %>% 
    data.frame() 
  
  cell_types <- rownames(bp_fractions) # extract cell types 
  rownames(bp_fractions) <- tolower(rownames(bp_fractions)) # lowercase for consistency
  
  bp_results[['fractions']] <- SummarizedExperiment(assays = list(cell_fractions = bp_fractions))
  
  # profiles ----
  if(is.na(key_cell) != TRUE){
    
    # retrieve estimated profiles for key cell (psi_mal)
    psi_mal <- bp.res@reference.update@psi_mal %>% 
      t() %>% 
      data.frame()
    
    # retrieve estimated profiles for the rest of cells (psi_env)
    psi_env <- bp.res@reference.update@psi_env %>% 
      t() %>% 
      data.frame()
    
    # join all profiles (group mode: 1 profile per cell type)
    psi_all <- cbind(psi_mal_avg = rowMeans(psi_mal), psi_env)
    colnames(psi_all) <- c(key_cell, colnames(psi_env))
    
    # add to list of results
    bp_results[['psi_mal']] <- SummarizedExperiment(assays = list(counts = psi_mal))
    bp_results[['geps']] <- SummarizedExperiment(assays = list(counts = psi_all))
    
  } else{
    psi_all <- bp.res@reference.update@phi %>% 
      t() %>% 
      data.frame()
    
    bp_results[['geps']] <- SummarizedExperiment(assays = list(counts = psi_all))
  }
  
  # arrange in SE and return ----
  bp_results <- MultiAssayExperiment(bp_results,
                                     metadata = list(tool = "BayesPrism",
                                                     key = key_cell))
  return(bp_results)
}
