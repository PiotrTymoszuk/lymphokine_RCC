# This script contains project globals

# libraries ----

  library(plyr)
  library(tidyverse)
  library(stringi)
  library(trafo)

# data container ------

  globals <- list()

# graphics -----

  globals$common_text <- element_text(size = 8,
                                      face = 'plain',
                                      color = 'black')

  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = 'mm')

  globals$common_theme <- theme_classic() + theme(axis.text = globals$common_text,
                                                  axis.title = globals$common_text,
                                                  plot.title = element_text(size = 8,
                                                                            face = 'bold'),
                                                  plot.subtitle = globals$common_text,
                                                  plot.tag = element_text(size = 8,
                                                                          face = 'plain',
                                                                          color = 'black',
                                                                          hjust = 0,
                                                                          vjust = 1),
                                                  plot.tag.position = 'bottom',
                                                  legend.text = globals$common_text,
                                                  legend.title = globals$common_text,
                                                  strip.text = globals$common_text,
                                                  strip.background = element_rect(fill = 'gray95',
                                                                                  color = 'gray80'),
                                                  plot.margin = globals$common_margin,
                                                  panel.grid.major = element_line(color = 'gray90'))

# genes of interest ------
  
  ## CCL4 and CCL5 are missing from the RECA-EU study

  globals$cxcl_genes <- c('CXCL9', 
                          'CXCL10', 
                          'CXCL11', 
                          'CXCR3', 
                          'CXCL13', 
                          'CXCR5', 
                          'XCL1', 
                          'XCR1')

# clinical variables -----

  globals$clin_vars <- c('sex', 'age', 'race', 'neoadjuvant',
                        'laterality', 'tumor_grade', 'p_stage',
                        'pt_stage', 'pm_stage', 'pn_stage',
                        'MSKCC_risk', 'IMDC_risk',
                        'sarc_diff', 'rhab_diff',
                        'tumor_delta_vol', 'response_type',
                        'benefit')

  globals$surv_vars <- c('relapse', 'rfs_days',
                         'death', 'os_days',
                         'tumor_death')

  globals$var_lexicon <- c(globals$var_labs,
                        set_names(c('Sex', 'Age, years', 'Race',
                                    'Neoadjuvant treatment', 'Laterality',
                                    'Tumor grade', 'Pathological stage',
                                    'Tumor stage', 'Metastasis stage',
                                    'Node stage, pn',
                                    'MSKCC risk group',
                                    'IMDC risk group',
                                    'Sarcomatoid differentiation',
                                    'Rhabdoid differentiation',
                                    'Tumor shrinkage',
                                    'Therapy response',
                                    'Clinical benefit'),
                                  globals$clin_vars),
                        set_names(c('Relapse', 'RFS, days',
                                    'Death', 'OS, days', 'Tumor-related death'),
                                  globals$surv_vars),
                        set_names(c('Overall response',
                                    'biresponse')))
  
  globals$var_lexicon <- globals$var_lexicon %>% 
    compress(names_to = 'variable', 
             values_to = 'label')

# infiltration estimates, cell types and labels -----
  
  ## quantiseq

  globals$quantiseq_lexicon <-
    c('T cell CD4+ (non-regulatory)' = 'CD4+ T',
      'T cell regulatory (Tregs)' = 'Treg',
      'T cell CD8+' = 'CD8+ T',
      'NK cell' = 'NK',
      'B cell' = 'B',
      'Myeloid dendritic cell' = 'mDC',
      'Neutrophil' = 'Neutro',
      'Monocyte' = 'Mono',
      'Macrophage M1' = 'M1 TAM',
      'Macrophage M2' = 'M2 TAM',
      'uncharacterized cell' = 'rest')
  
  globals$quantiseq_lexicon <- globals$quantiseq_lexicon %>% 
    compress(names_to = 'variable', 
             values_to = 'label')
  
  ## xCell: only the revelant cells are kept
  
  globals$xcell_lexicon <- 
    c('T cell gamma delta' = 'T \u03B3/\u03B4', 
      'T cell CD4+ (non-regulatory)' = 'CD4+ T', 
      'T cell CD4+ naive' = 'CD4+ Tn', 
      'T cell CD4+ memory' = 'CD4+ Tm', 
      'T cell CD4+ effector memory' = 'CD4+ Tem', 
      'T cell CD4+ central memory' = 'CD4+ Tcm', 
      'T cell CD4+ Th1' = 'CD4+ Th1', 
      'T cell CD4+ Th2' = 'CD4+ Th2', 
      'T cell regulatory (Tregs)' = 'Treg', 
      'T cell CD8+' = 'CD8+ T', 
      'T cell CD8+ naive' = 'CD8+ Tn', 
      'T cell CD8+ effector memory' = 'CD8+ Tem', 
      'T cell CD8+ central memory' = 'CD8+ Tcm', 
      'T cell NK' = 'NKT', 
      'NK cell' = 'NK', 
      'B cell' = 'B', 
      'B cell naive' = 'Bn', 
      'B cell memory' = 'Bm', 
      'Class-switched memory B cell' = 'switch Bm', 
      'B cell plasma' = 'B plasma', 
      'Myeloid dendritic cell' = 'mDC',
      'Myeloid dendritic cell activated' = 'act mDC', 
      'Plasmacytoid dendritic cell' = 'pDC', 
      'Neutrophil' = 'Neutro', 
      'Eosinophil' = 'Eosino', 
      'Mast cell' = 'MC', 
      'Monocyte' = 'Mono', 
      'Macrophage' = 'TAM', 
      'Macrophage M1' = 'M1 TAM', 
      'Macrophage M2' = 'M2 TAM', 
      'Endothelial cell' = 'EC', 
      'Cancer associated fibroblast' = 'CAF')

  globals$xcell_lexicon <- globals$xcell_lexicon %>% 
    compress(names_to = 'variable', 
             values_to = 'label')

# study cohorts ------

  globals$cohort_lexicon <- c(tcga = 'TCGA KIRC',
                              emtab1980 = 'E-MTAB 1980',
                              gse73731 = 'GSE73731',
                              gse167093 = 'GSE167093',
                              reca = 'RECA-EU',
                              cm10 = 'CM 010',
                              cm25ev = 'CM 025 EVER',
                              cm25ni = 'CM 025 NIVO')
  
  globals$cohort_lexicon <- globals$cohort_lexicon %>% 
    compress(names_to = 'cohort', 
             values_to = 'label')
  
  ## cohorts with survival
  
  globals$os_cohorts <- 
    c('tcga', 'emtab1980', 'reca', 'cm10', 'cm25ev', 'cm25ni')
  
# Lymphokine gene clusters -------
  
  globals$cluster_colors <- c('chemox low' = 'steelblue3', 
                              'chemox high' = 'coral3')

# END -----