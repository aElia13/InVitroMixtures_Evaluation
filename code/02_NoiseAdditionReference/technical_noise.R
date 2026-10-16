# Elia Alonso
# 17-12-2025

# load libraries
library("dplyr")
library("ModelMetrics")
library("ggplot2")

# Addition of white noise to reference

# some gene g
# count in profile C is Cg
# with volume (v) dependent noise

# Cg' ~ Pois(v*Cg) / v

# noise level ~ 1/v

# define function ----
add_technical_noise <- function(C, v){

  # ensure counts is a matrix
  C <- as.matrix(C) 
  
  # create empty matrix to store noise-corrupted profiles in the same positions
  C_noise <- matrix(0, nrow=nrow(C), ncol=ncol(C), dimnames = dimnames(C))
  
  # noise-level ~ 1/v - to make it directly proportional to v
  v <- 1 / v 
  
  for (i in 1:ncol(C)){ # iterate over profiles in the matrix to add noise
    # substract profile
    Ci <- C[,i]
    
    # calculate lambda of Poisson distribution
    lambda <- Ci * v
    # compute noisy counts
    Ci_noise <- rpois(n = length(lambda), lambda = lambda) / v
    # store noisy counts in matrix
    C_noise[,i] <- Ci_noise
  }
  
  # return noise-corrupted profile
  return(C_noise)
}

# # define working directory ----
# WD <- "/home/idc/elia.alonso/projects/invitro_mixtures"
# 
# # load data ----
# reference <- read.table(file.path("data/input_dig_cyt/invitro_pure_cell_lines_2025-08-25.txt")) 
# 
# # substract scc040 profiles from reference - to add noise
# reference.scc040 <- reference %>% 
#   dplyr::select(starts_with("SCC040"))
# 
# print(hist(log(as.matrix(reference.scc040)+1),
#            main="no noise",
#            ylim = c(0, 5000)))
# 
# head(reference.scc040)
# 
# # define values of v ----
# noise_volume <- c(5, 25, 50, 100, 200, 300)
# 
# # add noise to scc040 profiles and replace in reference matrix ----
# noise_added_ref <- list()
# 
# for(v in noise_volume){
# 
#   # add noise 
#   reference_noise <- add_technical_noise(C=reference.scc040, v=v)
#   
#   # update reference
#   reference$SCC040_1 <- reference_noise[, "SCC040_1"]
#   reference$SCC040_3 <- reference_noise[, "SCC040_3"]
#   
#   # save ref mat with added noise
#  noise_added_ref[[as.character(v)]] <- reference
#  }
# 
# # export noisy data ----
# file_name <- file.path(WD, "data/input_dig_cyt", 
#                        "invitro_pure_cell_lines_tech_noise_scc040.Rdata")
# noise_added_ref.technical <- noise_added_ref
# save(noise_added_ref.technical, file = file_name)
