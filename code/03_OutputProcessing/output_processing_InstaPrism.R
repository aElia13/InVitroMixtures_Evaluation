# script to process output data from deconvolution with 
# - InstaPrism

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

process_output_instaprism <- function(ip.res){
  key_cell <- ip.res@key
  ip_results = list() # create empty list to store processed results
  
  # fractions ----
  ip_fractions <- ip.res@theta %>%
      data.frame()
  
  cell_types <- rownames(ip_fractions) # extract cell types 
  rownames(ip_fractions) <- tolower(rownames(ip_fractions)) # save in lowercase for consistency

  ip_results[['fractions']] <- SummarizedExperiment(assays = list(cell_fractions = ip_fractions))
  
  # profiles ----

  # retrieve estimated profiles for key cell (psi_mal)
  psi_mal <- ip.res@psi_mal %>% 
    data.frame()
  
  # retrieve estimated profiles for the rest of cells (psi_env)
  psi_env <- ip.res@psi_env %>% 
    data.frame()
  
  # join all profiles , group mode
  psi_all <- cbind(psi_mal_avg = rowMeans(psi_mal), psi_env)
  colnames(psi_all) <- c(key_cell, colnames(psi_env))
  
  # add to list of results
  ip_results[['psi_mal']] <- SummarizedExperiment(assays = list(counts = psi_mal))
  ip_results[['geps']] <- SummarizedExperiment(assays = list(counts = psi_all))
    
  # arrange in SE and return
  ip_results <- MultiAssayExperiment(ip_results,
                                     metadata = list(tool = "InstaPrism",
                                                     key = key_cell))
  return(ip_results)
}
