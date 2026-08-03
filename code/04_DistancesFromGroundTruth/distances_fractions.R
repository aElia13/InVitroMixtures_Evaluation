# functions to calculate distances between deconvoluted fractions and ground truth

# Elia Alonso
# 25-07-2025

# load libraries
library("dplyr")
library("tidyr")
library("tibble")

# set working directory
WORKING_DIR <- "/home/idc/elia.alonso/projects/invitro_mixtures/"

# load functions
source("distances_equations.R")

# implement function
calculate_dist_fractions <- function(gt, est, tool){
  # gt: df containing true proportions per sample (cell_type x sample)
  # est: df containing estimated proportions by the tool (cell_type x sample)
  # tool: name of the digital cytometry tool

  print("---------------------------------------------------------------------")
  print("... calculating distances between true and estimated proportions")
  
  # 0. empty list to store results to return
  to_return <- list() 
  
  # 1. harmonize rows and cols from both gt and estimated proportions ----
  gt <- gt[rownames(est), colnames(est)] # rows = cell types; cols = samples
  
  if(!(unique(rownames(gt) == rownames(est))==T)){ ## double-check cells
    stop("-> cells NOT in same order in both true and estimated profiles")
  }
  
  if(!(unique(colnames(gt) == colnames(est))==T)){ ## double-check samples
    stop("-> samples NOT in same order in both true and estimated profiles")
  }
    
  # 2. calculate distances and store in df ----
  distances_df <- data.frame(row_names = colnames(est),
                             mix_id = colnames(est),
                             feature = "proportions",
                             tool=tool) %>% 
    column_to_rownames("row_names")
  
  corr <- sapply(1:ncol(est), function(i) {
    dist_corr(est=as.numeric(est[, i]),
              gt=as.numeric(gt[, i]))
  })
  names(corr) <- colnames(est)

  rmse <- sapply(1:ncol(est), function(i) {
    dist_rmse(est=as.numeric(est[, i]),
              gt=as.numeric(gt[, i]))
  })
  names(rmse) <- colnames(est)

  mae <- sapply(1:ncol(est), function(i) {
    dist_mae(est=as.numeric(est[, i]),
             gt=as.numeric(gt[, i]))
  })
  names(mae) <- colnames(mae)

  mink <- sapply(1:ncol(est), function(i) {
    dist_mink(est=as.numeric(est[, i]),
              gt=as.numeric(gt[, i]))
  })
  names(mink) <- colnames(est)

  weig <- sapply(1:ncol(est), function(i) {
    dist_weig(est=as.numeric(est[, i]),
              gt=as.numeric(gt[, i]),
              between_fractions = TRUE)
  })
  names(weig) <- colnames(est)

  distances_df <- cbind(distances_df,
                        corr,
                        rmse,
                        mae,
                        mink,
                        weig)
  
  # 3. return results ----
  to_return[["distances"]] <- distances_df
  return(to_return)
}
