# Survival analysis
#
# 1) Univariate survival analysis (Cox) with splines of expression values 
# of the genes of interest
#
# 2) Survival analysis (univariate Cox) done separately for the early 
# and late stage tumors
#
# 3) Survival in the lymphokine expression strata defined 
# by the the largest difference in survival - for the entire datasets 
# and early/late stage tumors
#
# 4) Multi-parameter overall an relapse-free survival with regularized Cox 
# regression employing the genes of interest as explanatory variables
#
# 5) Multi-parameter modeling of overall and relapse-free survival with 
# clinical parameters (age, sex, grade. stage) and gene expression variables

# tools ------

  library(plyr)
  library(tidyverse)
  library(exda)
  library(trafo)
  library(soucer)
  library(rlang)
  library(ggrepel)
  library(ggtext)
  library(furrr)
  library(ggExtra)

  library(survival)
  library(coxExtensions)
  library(splines)
  library(lmqc)
  library(kmOptimizer)
  library(glmnet)
  library(caret)
  library(survminer)

  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# survival scripts -------
  
  insert_msg('Survival scripts')
  
  c('./survival scripts/univariable_os.R', 
    './survival scripts/univariable_rfs.R', 
    './survival scripts/stages.R', 
    './survival scripts/stage_cutpoints.R', 
    './survival scripts/multi_os.R', 
    './survival scripts/multi_rfs.R', 
    './survival scripts/clinical_os.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# END ------
  
  insert_tail()