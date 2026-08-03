# Script to impolement functions to plot
## differential gene expression analysis
## within dig cytometry evaluation 

# Eli Alonso
# 24-02-2026

# Load libraries
library("ggplot2")
library("dplyr")
library("tibble")

DGEA <- function(geps1, geps2,
                 stderr1 = NULL, stderr2 = NULL, 
                 padj.cutoff = 0.1){
  
  # geps1: expression matrix of cell type 1 (genes x samples)
  # geps2: expression matrix of cell type 2 (genes x samples)
  
  # stderr1: standard error of expression values from cell type 1 
  # stderr2: standard error of expression values from cell type 2 
  
  # if stderr are not provided, they will be estimated with fitted model
  if(is.null(stderr1) & is.null(stderr2)){ # BayesPrism and InstaPrism
    print("no stderr info -> estimated with sd ~ mean")
    
    stderr.model <- function(x) 0.114*x + 1.1229e-06 
    stderr1 <- stderr.model(geps1)
    stderr2 <- stderr.model(geps2)
  }
  
  # 1. make sure all genes are in the same order 
  genes <- rownames(geps1) # genes in geps1 as reference
  geps2 <- geps2[genes,] %>% as.matrix()
  
  stderr1 <- stderr1[genes,] %>% as.matrix()
  stderr2 <- stderr2[genes,] %>% as.matrix()
  
  # 2. init empty matrices to store results of each comparison
  vBetaZ_mat  <- matrix(NA, nrow = nrow(geps1), ncol = ncol(geps1))
  ZPs_mat  <- matrix(NA, nrow = nrow(geps1), ncol = ncol(geps1))
  Zqvals_mat  <- matrix(NA, nrow = nrow(geps1), ncol = ncol(geps1))
  
  # 3. run diff. expression analysis, from:
  #-> https://pubmed.ncbi.nlm.nih.gov/31960376/ 
  #-> DOI: 10.1007/978-1-0716-0301-7_7 
  
  for (i in 1:ncol(geps1)){ # iterate over samples
    ## 3.1. calculate values
    vBetaZ <- sapply(1:nrow(geps1), 
                     function(j) (geps1[j,i] - geps2[j,i]) / sqrt(stderr1[j,i]^2 + stderr2[j,i]^2)
    )
    ZPs <- 2*pnorm(-abs(vBetaZ))
    Zqvals <- p.adjust(ZPs, method="BH")
    
    ## 3.2. save in matrices
    vBetaZ_mat[, i] <- vBetaZ
    ZPs_mat[, i] <- ZPs
    Zqvals_mat[, i] <- Zqvals
  }
  
  # 4. results
  res <- data.frame(genes = genes,
                    vBetaZ = vBetaZ_mat,
                    ZPs = ZPs_mat,
                    Zqvals = Zqvals_mat
                    )
  
  ## 4.1. filter results according to padj cutoff
  res.filt <- res %>%
    na.omit() %>% # remove null values
    dplyr::filter(Zqvals <= padj.cutoff) # filter diff. exp. genes
  
  cutoffs <- list(padj = padj.cutoff)
  
  ## list of differentially expressed genes
  diff_exp <- res.filt$genes
  
  ## 4.2. list of ranked list of genes for further GSEA - filtered by pval
  rank_list <- res$vBetaZ
  names(rank_list) <- res$genes # name list with gene names
  rank_list <- sort(rank_list, decreasing = TRUE) # rank by value
  rank_list <- na.omit(rank_list) # remove null values
  
  # 5. return results
  return(list(res = res, # raw results
              res.filt = res.filt, # filtered results by padj cutoff
              diff_exp = diff_exp, # list of differential expressed genes
              rank_list = rank_list, # ranked list of genes for GSEA
              cutoffs = cutoffs, # p.adjusted cutoff used for filtering DEGs
              ) 
         )
}
