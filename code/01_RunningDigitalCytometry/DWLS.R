# script to run deconvolution with DWLS

# Elia Alonso
# 14-10-2025

# load libraries
library("DWLS")
library("dplyr")
library("tibble")

# define working directory
WD <- "/Users/elia.alonso/Documents/Code/invitro_mixtures/"
# WD <- "/home/idc/elia.alonso/projects/invitro_mixtures"

# load data ----
reference <- read.table(file.path(WD, "data/input_dig_cyt/invitro_pure_cell_lines_2025-08-25.txt")) 
mixtures <- read.table(file.path(WD, "data/input_dig_cyt/invitro_mixtures_bulk_data_2025-08-25.txt"))


### build and save signature matrix ----
pval_cutoff = 0.1 # for building signature matrix

dataSC <- reference %>%
  as.matrix()

# cell types labels:
labels <- c(rep("scc040", 2), rep("fibs", 2), rep("tcells", 2), rep("nk", 2))

res_path <- file.path(WD, "data/output_dig_cyt/DWLS/sig_mat_scc040_fibs_nk_tcells")
Signature <- buildSignatureMatrixMAST(scdata = dataSC,
                                      id = labels,
                                      path = res_path,
                                      diff.cutoff = 0.5, # default
                                      pval.cutoff = pval_cutoff) 

# save signature matrix
file_name <- file.path(res_path, "sig_mat_scc040_fibs_nk_tcells.csv")
write.csv(data.frame(Signature) %>% rownames_to_column("genes"), 
          file=file_name,
          row.names = F)

### run DWLS ----
DWLS_res <- list() # empty list to store results

mixtures.tmp <- mixtures %>% 
  data.frame() %>% 
  na.omit() # ensure that there are no missing entries

# iterate over samples and estimate
for (sample in colnames(mixtures.tmp)){ 
  print(sample)
  
  # extract counts of sample and assign gene symbols as names
  dataBulk <- mixtures.tmp[, sample]
  names(dataBulk) <- rownames(mixtures.tmp)
  
  # trimm signature and bulk data to contain same Diff. exp. genes
  tr <- trimData(Signature, dataBulk)
  S <- tr$sig # trimmed signature
  B <- tr$bulk # trimmed bulk
  
  # estimate using dampened weighted least squares
  solDWLS <- solveDampenedWLS(S, B)
  
  # store in list of results
  DWLS_res[[sample]] <- solDWLS
}

### export results ----

# transform into data frame (rows = cell types, columns = samples)
DWLS_res.df <- as.data.frame(DWLS_res)

# save
file_name <- file.path(WD, "data/output_dig_cyt/DWLS/DWLS_output_fractions_Ia.csv")
write.csv(DWLS_res.df, 
          file=file_name,
          row.names = T)
