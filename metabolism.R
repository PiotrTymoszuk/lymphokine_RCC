# Differences in metabolism based on differential gene expression
# between the lymphokine clusters

# Tools -------

  library(plyr)
  library(tidyverse)
  library(exda)
  library(microViz)
  library(soucer)
  library(ggrepel)
  library(stringi)
  library(furrr)
  library(BiGGR)
  library(biggrExtra)
  library(rstatix)
  library(ggtext)

  explore <- exda::explore
  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)

# scripts ------

  insert_msg('Metabolism scripts')
  
  ## construction of metabolism models - resorting to the cache
  
  if(file.exists('./cache/metabolism.RData')) {
    
    load('./cache/metabolism.RData')
    
  } else {
    
    source_all('./metabolism scripts/testing.R', 
               message = TRUE, crash = TRUE)
    
  } 

  c('./metabolism scripts/plots.R', 
    './metabolism scripts/detail_plots.R', 
    './metabolism scripts/gene_plots.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
