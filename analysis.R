# The following analyses are done
#
# 1) Characteristic of the study cohorts
#
# 2) Screening of correlations of all chemokine genes with predicted CD8+ T cell 
# content (Quantiseq) and of regulation of chemokine genes in tumor vs normal 
# kidney. This step is done only in the TCGA cohort, our primary screening 
# collective
#
# 3) Normal - RCC comparison of the CXCL9/10/11 and CXCR3 expression, TCGA
# and GSE167093
#
# 4) Co-expression of CXCL9, CXCL10, CXCL11 and CXCR3 in all study cohorts
#
# 5) Expression of the genes of interest in tumor stages and demographic strata
#
# 6) Correlation of the genes of interest with QuanTIseq infiltration estimates
#
# 7) Correlation of the genes of interest with xCell infiltration estimates
#
# 8) UMAP of by the xCell infiltration estimates and analysis 
# of overlap between the gene expression and infiltration
#
# 9) Multi-dimensional scaling to assess association of the genes of interest
# with the xCell infiltration estimates


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
  library(dbscan)
  library(igraph)
  library(network)
  library(ggnetwork)

  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select

  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# analysis scripts ------
  
  insert_msg('Analysis scripts')
  
  c('./analysis scripts/cohorts.R', 
    './analysis scripts/screening.R', 
    './analysis scripts/norm_tumor.R', 
    './analysis scripts/coexpression.R', 
    './analysis scripts/clinics.R', 
    './analysis scripts/quantiseq.R', 
    './analysis scripts/xcell.R', 
    './analysis scripts/clustering.R', 
    './analysis scripts/mds.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# END ------
  
  insert_tail()