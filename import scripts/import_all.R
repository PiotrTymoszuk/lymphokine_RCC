# This script imports clinical, expression, gene annotation
# infiltration (quanTIseq, xCell) and signature data (Reactome, Recon)

# toolbox ----

  library(plyr)
  library(tidyverse)
  library(stringi)
  library(soucer)
  library(xena)
  library(readxl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(furrr)
  library(GEOquery)
  library(immunedeconv)
  library(rlang)
  library(trafo)
  
  library(biggrExtra)
  library(gseaTools)

  reduce <- purrr::reduce
  select <- dplyr::select

  insert_head()

# data containers -----

  tcga_expr <- list()
  tcga_clinics <- list()

  somut <- list() ## mutations

  tcga <- list() ## for the final TCGA data container

  ## containers for validation data sets

  emtab1980 <- list()
  gse73731 <- list()
  gse167093 <- list()
  cm10 <- list() ## Checkmate 10
  cm25ev <- list() ## checkmate 25 evrolimus
  cm25ni <- list() ## checkmate 25 nivolumab
  reca <- list() ## EU-RECA
  
  ## containers for deconvolution results
  ## Reactome and Recon2 signatures
  
  deconv <- list() ## infiltration estimates
  sig <- list() ## Recon and Reactome signatures

# globals setup -----

  insert_msg('Globals setup')

  source_all('./tools/project_globals.R',
             message = TRUE, crash = TRUE)

# clearing the clinical and expression data sets ------

  insert_msg('Clearning the clinical and expression data')

  c('./import scripts/data_clearing_tcga_clinics.R',
    './import scripts/data_clearing_tcga_expression.R',
    './import scripts/emtab1980.R',
    './import scripts/gse73731.R',
    './import scripts/gse167093.R',
    './import scripts/checkmate.R',
    './import scripts/reca.R') %>%
    source_all(source, message = TRUE, crash = TRUE)

# merging the expression data with the essential set of clinical information -----

  insert_msg('Appending the expression data with the essential clinical information')

  tcga[c('clinical',
         'drugs',
         'radiation')] <- tcga_clinics[c('clinical',
                                         'drugs',
                                         'radiation')]

  tcga[c('expression',
         'annotation')] <- tcga_expr[c('exprs',
                                       'annotation')]

  tcga$expression <- right_join(tcga$clinical,
                                tcga$expression,
                                by = 'patient_id')

# Duplicate handling: samples with duplicated analyte ID are removed ------

  insert_msg('Duplicate handling, removing duplicated analyte id')

  ## duplicated analyte id

  tcga$expression <- tcga$expression %>%
    dlply(.(tissue_type)) %>%
    map_dfr(filter, !duplicated(analyte_id)) %>%
    as_tibble

# XENA mutations ------

  insert_msg('XENA mutations')

  source_all('./import scripts/xena.R',
             crash = TRUE, message = TRUE)

# Infiltration estimates ------

  insert_msg('Infiltration estimates')

  ## working primarily with cache
  ## the TCGA and GSE167093 infiltration datasets contain
  ## the tumor samples only

  if(!file.exists('./data/infiltration/deconv.RData')) {

    source_all('./import scripts/infiltration.R',
               message = TRUE, crash = TRUE)

  } else {

    load('./data/infiltration/deconv.RData')

  }

  deconv <- deconv %>%
    map(~map(.x, 
             column_to_rownames, 
             'cell_type') %>% 
          map(t) %>%
          map(as.data.frame) %>%
          map(rownames_to_column, 'patient_id') %>%
          map(as_tibble))
  
  ## appending the study data lists
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0(i, '$quantiseq <- deconv$quantiseq_results$', i)))
    eval(parse_expr(paste0(i, '$xcell <- deconv$xcell_results$', i)))
    
  }
  
# Recon and Reactome signatures ----------
  
  insert_msg('Recon and Reactome signatures')
  
  ## external script for calculation of the signatures
  ## or cache
  
  if(file.exists('./data/signatures/signatures.RData')) {
    
    load('./data/signatures/signatures.RData')
    
  } else {
    
    source_all('./import scripts/signatures.R', 
               message = TRUE, crash = TRUE)
    
  }
  
  ## appending the study data lists
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0(i, '$reactome <- sig$reactome_scores$', i)))
    eval(parse_expr(paste0(i, '$recon <- sig$recon_scores$', i)))
    
  }
  
  sig <- sig[c('reactome_lexicon', 
               'recon_lexicon')]

# saving the workspace -------

  insert_msg('Saving the workspace')

  save(tcga, 
       emtab1980, gse73731, gse167093, 
       reca, cm10, cm25ev, cm25ni,
       somut, sig, 
       file = './data/processed_data.RData')

# END ----

  rm(tcga_clinics,
     tcga_expr,
     deconv,
     i)

  insert_tail()
