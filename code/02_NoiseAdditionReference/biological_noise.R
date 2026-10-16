# Elia Alonso
# 24-11-2025

# load libraries
library("dplyr")
library("ModelMetrics")
library("tidyr")
library("ggplot2")

# Addition of biological noise to reference

# to mimic cell cycle states, interactions, influence of TME on transc. prof...

# define function ----
add_biological_noise <- function(base_cell, add_cell, v = 0.1){
  # base_cell = profile of cells to add noise on
  # add_cell = profile of cells that will be added to base_cell
  # v = volume [0, 1] in which add_cell is added to base_cell
  
  # compute noisy counts
  C_noise <- (1 - v)*base_cell + v*add_cell
  
  # return noisy profile
  return(C_noise)
}

# define working directory ----
WD <- "/home/idc/elia.alonso/projects/invitro_mixtures"

# load data ----
reference <- read.table(file.path(WD, "data/input_dig_cyt/invitro_pure_cell_lines_2025-08-25.txt")) 

# substract scc040 profiles from reference - to add noise on
reference.scc040 <- reference %>% 
  dplyr::select(starts_with("SCC040"))

# substract okf6 profiles from reference - to mix with SCC040
reference.okf6 <- reference %>% 
  dplyr::select(starts_with("OKF6"))


# define noise levels ----
noise_volume <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1)

# add noise to scc040 profiles and replace in reference matrix ----
noise_added_ref <- list()

for(v in noise_volume){
  reference_noise <- add_biological_noise(base_cell = reference.scc040, 
                                          add_cell = reference.okf6,
                                          v=v)
  
  # update reference with corrupted profiles of SCC040, mixed with OKF6
  reference$SCC040_1 <- reference_noise[, "SCC040_1"]
  reference$SCC040_3 <- reference_noise[, "SCC040_3"]
  
  # save ref mat with added noise
  noise_added_ref[[as.character(v)]] <- reference
}


# export ----
file_name <- file.path(WD, "data/input_dig_cyt",
                       "invitro_pure_cell_lines_biol_noise_scc040.Rdata")
noise_added_ref.biological <- noise_added_ref
save(noise_added_ref.biological, file = file_name)

