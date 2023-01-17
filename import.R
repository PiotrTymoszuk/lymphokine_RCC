# This script imports TCGA RNA Seq and clinical data, works with cache.

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(soucer)
  library(rlang)

  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()

  source_all('./tools/project_globals.R',
             message = TRUE, crash = TRUE)

# loading the processed data ------

  insert_msg('Loading the processed data')

  if(file.exists('./data/processed_data.RData')) {

    load('./data/processed_data.RData')

  } else {

    source_all('./import scripts/import_all.R')

  }

# END ----

  insert_tail()
