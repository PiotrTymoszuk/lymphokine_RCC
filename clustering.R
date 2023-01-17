# Semi-supervised clustering of the tumor samples by the genes of interest
#
# The clusters are trained in the TCGA cohort and predicted for the remaining 
# collectives with a kNN classifier. The following analyses are done for 
# the lymphokine clusters
#
# 1) comparison of predicted immune infiltration levels (xCell and QuanTIseq)
#
# 2) comparison of overall and relapse-free survival 
# (uni-variable/log-rank and multi-variate Ridge Cox)
#
# 3) comparison of clinical features (demographics, grade and stages)
#
# 4) Differences in expression of metabolic gene signatures (Recon2 subsystems)
# and Reactome biological process gene signatures
#
# 5) Differences in mutation burden and gene mutation frequency between 
# the clusters

# tools -------

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
  library(somKernels)

  library(survival)
  library(survminer)
  library(glmnet)
  library(coxExtensions)

  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# clustering scripts ------
  
  insert_msg('Clustering scripts')
  
  ## finding the optimal gene clustering algorithm
  ## working with cache - time costly
  
  if(file.exists('./cache/gene_clusters.RData')) {
    
    load('./cache/gene_clusters.RData')
    
  } else {
    
    source_all('./clustering scripts/cluster_development.R', 
               message = TRUE, crash = TRUE)
    
  }

  c('./clustering scripts/clustering.R', 
    './clustering scripts/xcell.R', 
    './clustering scripts/quantiseq.R', 
    './clustering scripts/survival.R', 
    './clustering scripts/multi_survival.R', 
    './clustering scripts/clinics.R', 
    './clustering scripts/recon.R', 
    './clustering scripts/reactome.R', 
    './clustering scripts/mutations.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END -------
  
  insert_tail()