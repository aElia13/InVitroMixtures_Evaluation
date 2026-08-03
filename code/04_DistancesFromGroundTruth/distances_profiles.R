# functions to calculate distances between deconvoluted profiles and ground truth

# Elia Alonso
# 20-02-2026

# load libraries
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("ggExtra")
library("viridis")
library("ggridges")
library("gridExtra")
library("grid")

# set working directory
WORKING_DIR <- "/home/idc/elia.alonso/projects/invitro_mixtures/"

# load functions ----
source("distances_equations.R")

# load data ----
## pre-defined bins to split profiles and calculate distances bin-wise
## available in Supplementary Table 2.

bins <- read_table(file = file.path(WORKING_DIR, "data/bins_exp_genes.csv"))

# implement function ----
calculate_dist_profiles_group_mode <- function(gt, est, bins, tool){
  # gt: SE containing true profiles per pure sample (genes x sample)
  # est: df containing estimated profiles by the tool (genes x cell_types)
  # bins: list where names are thresholds names and values thresholds values
  # tool: name of the digital cytometry tool

  # 0. create empty objects to store data ----
  absolute.deviation <- data.frame() # empty df to store abs. dev. bin-gene-wise
  distances <- data.frame() # empty df to store distances bin-wise

  # 0.1. check that the tool name is correct ----
  supported_tools <- c("CIBERSORTx", "BayesPrism", "InstaPrism")
  if(!(tool %in% supported_tools)){
    print("please specify the name of one of the supported tools:")
    print(supported_tools)
    stop("-> tool name not present in supported tools")
  } 
  
  # 0.2. function to normalize to 1 profiles before distance calculation ----
  norm_to_1 <- function(x){ # L1 normalization
    if (sum(x, na.rm = TRUE) == 0) { # if there is no signal detected --> 0
      norm_x <- x / 1
    } else { # if there is signal detected --> norm to sum up to 1
      norm_x <- x / sum(x, na.rm = TRUE)
    }
    return(norm_x) 
  }
  
  # 0.3. retrieve profiles of pure samples (ground truth) ----
  gt_samples_cells <- colData(gt) %>% 
    data.frame() %>% 
    dplyr::select(c(sample_id, cell_id))
  
  gt_samples_cells <- gt_samples_cells %>% # create cell_type to match bins info
    dplyr::mutate(cell_type = gsub("_.*", "", sample_id)) %>% 
    dplyr::mutate(cell_type = ifelse(cell_type=="Fibs", "Fib", cell_type),
                  cell_type = ifelse(cell_type=="Tcells", "TC", cell_type))
  gt <- assay(gt)
  
  for (ct in colnames(est)){ # iterate over cell types
    
    # 1. retrieve ground truth and estimated profiles from this cell type ----
    ## ground truth
    ct.samples <- gt_samples_cells %>% 
      dplyr::filter(cell_id==ct) %>% 
      dplyr::pull(sample_id)
    gt.ct <- gt[, ct.samples]
    
    ## estimated by tool
    est.ct <- est %>% 
      dplyr::select(all_of(ct)) %>% 
      data.frame()
    
    # 2. filter 1 and NA from Cx estimated profile ----
    if (tool == "CIBERSORTx"){ # remove 1 and NA of estimated profile of each cell type
      # The "1" values in the expression matrix txt files are genes with insufficient
      # evidence of expression (these genes are either not expressed or have inadequate
      #                         statistical power to be imputed).
      # The NA values are genes that have inadequate statistical power to be imputed.
      
      est.ct <- est.ct %>%
        dplyr::filter(.data[[ct]] != 1) %>% # remove genes with 1
        na.omit() # remove null values
      }
    
    # 3. subset genes in common between true and estimated profiles ----
    common_genes <- intersect(rownames(gt.ct), rownames(est.ct))
    gt.ct <- gt.ct[common_genes,] # subset from ground truth
    est.ct <- est.ct[common_genes, , drop=F] # subset from estimated profile
    
    # 4. scale profiles to sum 1 (L1 normalization) ----
    gt.ct_norm <- as.data.frame(apply(gt.ct, 2, norm_to_1))
    
    if (tool %in% c("BayesPrism", "InstaPrism")){ # already normalized
      est.ct_norm <- est.ct
    } else{ # normalize profiles with function defined above
      est.ct_norm <- as.data.frame(apply(est.ct, 2, norm_to_1))
    }
    
    if((unique(rownames(gt.ct_norm)==rownames(est.ct_norm))==F)){ # check genes
      stop("-> genes NOT in same order in both true and estimated profiles")
    }
    
    # 5. bins of expression ----
    ## retrieve cell type id from the ground truth profile to retrieve bins
    cell_type <- gt_samples_cells %>% 
      dplyr::filter(cell_id == ct) %>% 
      pull(cell_type) %>% 
      unique()
    
    # retrieve bin info for the cell type ground truth
    ct.bin <- bins %>% 
      dplyr::filter(genes %in% common_genes) %>% # subset genes present in profs
      dplyr::select(c("genes", starts_with(cell_type)&ends_with("bin"))) 
    ct.bin <- split(ct.bin[[1]], factor(ct.bin[[2]]))
    
    for (bin in names(ct.bin)){ # iterate over bins
      print(bin)
      genes_bin <- ct.bin[[bin]]
      n = length(genes_bin)

      if(n==0){ # no genes in this bin --> all results = NA
        df <- data.frame(
          corr = NA,
          rmse = NA,
          mae = NA,
          cell_type = ct,
          expression_level = bin,
          feature = "profiles",
          tool = tool,
          n_genes = n
        )
        distances <- rbind(distances, df)
        
      } else{ # subset norm profiles and calculate distances
        gt.ct_norm.bin <- gt.ct_norm[genes_bin, ct.samples] 
        est.ct_norm.bin <- est.ct_norm[genes_bin, ct] 
        
        # empty df to store results from both replicates
        tmp <- list(abs = data.frame(row.names = genes_bin),
                    cor = c(),
                    rms = c(),
                    mae = c())
        
        for(i in 1:ncol(gt.ct_norm.bin)){ # iterate over replicates of pure profs
          gt.ct_norm.bin.i <- gt.ct_norm.bin[,i]
          
          # calculate abs dev (gene-wise)
          tmp$abs[,i] <- abs(gt.ct_norm.bin.i - est.ct_norm.bin)
          
          # calculate distances: 1-p.corr, rmse, mae (bin-wise)
          tmp$cor <- c(tmp$cor, dist_corr(gt = gt.ct_norm.bin.i, 
                                          est = est.ct_norm.bin))
          tmp$rms <- c(tmp$rms, dist_rmse(gt = gt.ct_norm.bin.i, 
                                          est = est.ct_norm.bin))
          tmp$mae <- c(tmp$mae, dist_mae(gt = gt.ct_norm.bin.i, 
                                         est = est.ct_norm.bin))
        }
        
        # merge absolute deviations
        dev_results <- data.frame(
          genes = rownames(tmp$abs),
          abs_dev = rowMeans(tmp$abs)
        )
        dev_results$cell_type <- ct
        dev_results$expression_level <- bin
        dev_results$tool <- tool
        
        absolute.deviation <- rbind(absolute.deviation, dev_results)
        
        # merge distances
        dist_results <- data.frame(
          row.names = c("corr", "rmse", "mae"),
          values = c(mean(tmp$cor),
                      mean(tmp$rms),
                      mean(tmp$mae)
                      )
        )
      
        dist_results <- as.data.frame(t(dist_results))
        dist_results$cell_type = ct
        dist_results$expression_level = bin
        dist_results$feature = "profiles"
        dist_results$tool = tool
        dist_results$n_genes = n

        distances <- rbind(distances, dist_results)
      }
      
    }
  }
  
  # 6. return ----
  distances <- distances %>%
    dplyr::select(c(cell_type, feature, tool, 
                    expression_level, n_genes,
                    corr, rmse, mae))
  
  return(list("distances" = distances,
              "abs.deviation" = absolute.deviation)
         )
}
