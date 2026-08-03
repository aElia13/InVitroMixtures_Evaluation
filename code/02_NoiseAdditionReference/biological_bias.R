# Elia Alonso
# 24-11-2025

# load libraries
library("dplyr")
library("ModelMetrics")
library("tidyr")
library("ggplot2")

# Addition of biological bias to reference

# to mimic cell cycle states, interactions, influence of TME on transc. prof...

# define function ----
add_biological_bias <- function(base_cell, add_cell, v = 0.1){
  # base_cell = profile of cells to add bias on
  # add_cell = profile of cells that will be added to base_cell
  # v = volume [0, 1] in which add_cell is added to base_cell
  
  # compute noisy counts
  C_bias <- (1 - v)*base_cell + v*add_cell
  
  # return noisy profile
  return(C_bias)
}


# # load data ----
# reference <- read.table(file.path("data/input_dig_cyt/invitro_pure_cell_lines_2025-08-25.txt")) 
# 
# # substract scc040 profiles from reference - to add bias on
# reference.scc040 <- reference %>% 
#   dplyr::select(starts_with("SCC040"))
# 
# # substract okf6 profiles from reference - to mix with SCC040
# reference.okf6 <- reference %>% 
#   dplyr::select(starts_with("OKF6"))
# 
# 
# # define bias levels ----
# bias_volume <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1)
# 
# # add increasing levels of bias to scc040 profiles ----
# bias_added_ref <- list()
# 
# for(v in bias_volume){
#   reference_bias <- add_biological_bias(base_cell = reference.scc040, 
#                                           add_cell = reference.okf6,
#                                           v=v)
#   
#   # update reference with corrupted profiles of SCC040, mixed with OKF6
#   reference$SCC040_1 <- reference_bias[, "SCC040_1"]
#   reference$SCC040_3 <- reference_bias[, "SCC040_3"]
#   
#   # save ref mat with added bias
#   bias_added_ref[[as.character(v)]] <- reference
# }
