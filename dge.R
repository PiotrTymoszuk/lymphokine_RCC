# Differential gene expression between the lymphokine clusters
#
# 1) Differentially regulated genes (all and exhaustion-related) are identified
# by two-tailed FDR-corrected T test and at least 1.5-fold regulation between 
# the clusters
#
# 2) GO enrichment for the differentially regulated genes.
#
# 3) Differential gene expression for the selected functional groups of genes

# tools --------

  library(plyr)
  library(tidyverse)
  library(exda)
  library(trafo)
  library(soucer)
  library(rlang)
  library(microViz)
  library(ggrepel)
  library(ggtext)
  library(furrr)
  library(clustTools)
  library(readxl)

  library(limma)
  library(SPIA)

  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# differential gene expression scripts -------
  
  insert_msg('Differential gene expression scripts')
  
  ## serial testing is time costly and done via cache
  
  if(file.exists('./cache/dge.RData')) {
    
    load('./cache/dge.RData')
    
  } else {
    
    source_all('./dge scripts/expression.R', 
               message = TRUE, crash = TRUE)
    
  }
  
  c('./dge scripts/plots.R', 
    './dge scripts/go.R', 
    './dge scripts/exhaustion.R', 
    './dge scripts/interferon.R', 
    './dge scripts/il9_21.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
  ## SPIA results from the cache (extremely computation intensive)

  if(file.exists('./cache/spia.RData')) {
    
    
    load('./cache/spia.RData')
    
  } else {
    
    source_all('./dge scripts/spia.R', 
               message = TRUE, crash = TRUE)
    
  }
    
# END --------
  
  insert_tail()