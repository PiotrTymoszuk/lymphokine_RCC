# Figures, supplementary material and manuscript

# tools -------

  library(plyr)
  library(tidyverse)
  library(soucer)
  library(stringi)
  library(furrr)

  library(exda)
  library(clustTools)

  library(cowplot)
  library(writexl)
  library(figur)
  library(rmarkdown)
  library(knitr)
  library(bookdown)
  library(flextable)
  
  explore <- exda::explore
  extract <- clustTools::extract
  nobs <- clustTools::nobs
  reduce <- purrr::reduce
  select <- dplyr::select
  
  insert_head()
  
  c('./tools/project_globals.R', 
    './tools/project_tools.R') %>% 
    source_all(message = TRUE, crash = TRUE)
  
# paper scripts -------
  
  insert_msg('Paper scripts')
  
  c('./paper scripts/tables.R', 
    './paper scripts/supplement.R', 
    './paper scripts/figures.R', 
    './paper scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)
  
# END -------
  
  insert_tail()