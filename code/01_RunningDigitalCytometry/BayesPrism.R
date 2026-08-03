#
# Script to run BayesPrism
#

# elia.alonso@helmholtz-munich.de
# 22-09-2025

# This script:
# - Loads input data (must be a Prism object)
# - Checks that input data is a Prism object
# - Runs BayesPrism on input prism object 
# - Exports results in .RDS file

#----Load libraries----
library("BayesPrism")

#----Variables----
WORKING_DIR = "/home/idc/elia.alonso/projects/invitro_mixtures" # ! Note: adjust 
setwd(WORKING_DIR)

INPUT_PATH = "./data/input_dig_cyt/Bayes_Insta_Prism/"
files <- list.files(INPUT_PATH, pattern = "\\.Rdata$", full.names = TRUE)

OUTPUT_DATA = "./data/output_dig_cyt/BayesPrism/BayesPrism_results.Rdata"

# create empty list to store BayesPrism results
bp_results <- list()

for (INPUT_DATA in files){
  res <- list() # list to store deconv results and time of execution
  
  #----Load input data----
  print(paste("...loading", INPUT_DATA))
  filename <- tools::file_path_sans_ext(basename(INPUT_DATA))

  # check that input file exists
  if (!file.exists(INPUT_DATA)) {
    stop("Error: File does not exist.")
  }
  
  # check that input file contains information
  if (is.na(file.info(INPUT_DATA)$size) || file.info(INPUT_DATA)$size == 0) {
    stop("Error: File is empty or corrupted.")
  }
  
  # load input data
  loaded_object <- load(INPUT_DATA) # load data object
  
  assign("input_prism", get(loaded_object[1])) # assign fixed variable name for consistency
  rm(list = loaded_object) # remove duplicated data for memory
  rm(loaded_object)
  
  # Check that data format is Prism (required by BayesPrism)
  if (!inherits(input_prism, "prism")) {
    stop("Error: The input object is not of class 'prism'.")
  } else {
    print("... input data loaded successfully.")
  }
  
  #----Run analysis----
  print("... running BayesPrism")
  
  bp.res <- run.prism( # run BayesPrism
    prism = input_prism,
    n.cores = 50
    )
  
  #----Store results----
  res[["bp_res"]] = bp.res
  bp_results[[filename]] = res
}

#----Export results----
print(paste("... saving results to file:", OUTPUT_DATA))
save(bp_results, file = OUTPUT_DATA) # save BayesPrism results
