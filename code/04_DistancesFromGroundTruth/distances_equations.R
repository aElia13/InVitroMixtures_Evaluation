# functions to calculate distances 
# between digital cytometry output and ground truth

# Elia Alonso
# 25-07-2025

# load libraries
library("stats")
library("ModelMetrics")

# 1 - p.correlation
dist_corr <- function(gt, est) {
  d_cor <- 1 - stats::cor(gt, est, method = "pearson")
  return(d_cor)
}

# root mean squared error (RMSE)
dist_rmse <- function(gt, est) {
  d_rmse <- ModelMetrics::rmse(actual=gt, predicted=est)
  return(d_rmse)
}

# mean absolute error (MAE)
dist_mae <- function(gt, est){
  d_mae <- ModelMetrics::mae(actual=gt, predicted=est)
}

# minkowski p=0.5
dist_mink <- function(gt, est, p=0.5) {
  d_mink <- sum(abs(est - gt)^p)^(1/p)
  return(d_mink)
}

# weighted distances
dist_weig <- function(gt, est, between_fractions){
  
  if (between_fractions == T){ # distances between fractions
    epsilon <- 0.05
  } else{ # distances between profiles
    epsilon <- log(10, base = 2)}
  
  w = 1/(gt + epsilon)
  d_wgt <- sum(w*abs(est - gt))
  return(d_wgt)
}
