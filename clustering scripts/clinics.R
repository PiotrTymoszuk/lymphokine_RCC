# Comparison of clinincal features between the lymphokine clusters

  insert_head()
  
# container -----
  
  clust_clinic <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## variables 
  
  clust_clinic$variables <- globals$clin_vars
  
  ## analysis tables with the clinical variables
  
  clust_clinic$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('clust_clinic$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  clust_clinic$analysis_tbl[c('tcga', 'gse167093')] <- 
    clust_clinic$analysis_tbl[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  clust_clinic$analysis_tbl <- clust_clinic$analysis_tbl %>% 
    map(~select(.x,
                patient_id, 
                any_of(clust_clinic$variables)))
  
  ## appending them with the clustering information
  
  clust_clinic$analysis_tbl <- gene_clust$clust_obj %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(., clust_clinic$analysis_tbl, 
         left_join, by = 'patient_id')
  
  ## manual re-coding of the grade and stage information: any TX, MX and NX
  ## are re-coded to NAs. Tumor stages are simplified (T1 - T4)
  
  clust_clinic$analysis_tbl[c('tcga', 'emtab1980', 'reca')] <- 
    clust_clinic$analysis_tbl[c('tcga', 'emtab1980', 'reca')] %>% 
    map(mutate, 
        tumor_grade = stri_extract(tumor_grade, regex = '\\d{1}'), 
        tumor_grade = paste0('G', tumor_grade), 
        tumor_grade = factor(tumor_grade, c('G1', 'G2', 'G3', 'G4')), 
        pt_stage = stri_extract(pt_stage, regex = 'T\\d{1}'), 
        pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')), 
        pm_stage = factor(pm_stage, c('M0', 'M1')), 
        pn_stage = factor(pn_stage, c('N0', 'N1', 'N2')))
  
  clust_clinic$analysis_tbl[c('gse73731', 'gse167093')] <- 
    clust_clinic$analysis_tbl[c('gse73731', 'gse167093')] %>% 
    map(mutate, 
        tumor_grade = stri_extract(tumor_grade, regex = '\\d{1}'), 
        tumor_grade = paste0('G', tumor_grade), 
        tumor_grade = factor(tumor_grade, c('G1', 'G2', 'G3', 'G4')), 
        pt_stage = stri_extract(pt_stage, regex = 'T\\d{1}'), 
        pt_stage = factor(pt_stage, c('T1', 'T2', 'T3', 'T4')))
  
  ## manual re-coding of the therapy response and tumor shrinkage information
  
  clust_clinic$analysis_tbl[c('cm10', 'cm25ev', 'cm25ni')] <-
    clust_clinic$analysis_tbl[c('cm10', 'cm25ev', 'cm25ni')] %>%
    map(mutate,
        tumor_delta_vol = ifelse(tumor_delta_vol < 0,
                                 'yes', 'no'),
        tumor_delta_vol = factor(tumor_delta_vol, 
                                 c('no', 'yes')),
        biresponse = car::recode(as.character(response_type),
                                 "'CR' = 'CR/PR';
                                 'PR' = 'CR/PR';
                                 'SD' = 'SD/PD';
                                 'PD' = 'SD/PD';
                                 'CRPR' = 'CR/PR'"),
        biresponse = factor(biresponse, 
                            c('SD/PD', 'CR/PR')))
  
  ## complete cases
  
  clust_clinic$analysis_tbl <- clust_clinic$analysis_tbl %>% 
    map(filter, !is.na(clust_id))
  
  ## independent clinical variables
  ## there's only one arm in the CheckMate 010 cohort,
  ## an rhabdoid differentiation
  ## the variable is hence removed
  
  clust_clinic$variables <- clust_clinic$analysis_tbl %>%
    map(names) %>%
    map(~.x[.x %in% c(globals$clin_vars, 'biresponse')])
  
  clust_clinic$variables$cm10 <-
    clust_clinic$variables$cm10[!clust_clinic$variables$cm10 %in% c('rhab_diff')]
  
  ## effect size statistics, dependent on the variable type
  ## plot_type and y label
  
  clust_clinic$types <-
    map2(clust_clinic$analysis_tbl,
         clust_clinic$variables,
         ~map_lgl(.x[.y], is.numeric)) %>% 
    map(~ifelse(.x, 'wilcoxon_r', 'cramer_v'))
  
  clust_clinic$plot_types <-
    map2(clust_clinic$analysis_tbl,
         clust_clinic$variables,
         ~map_lgl(.x[.y], is.numeric)) %>% 
    map(~ifelse(.x, 'violin', 'stack'))
  
  clust_clinic$y_labels <-
    map2(clust_clinic$analysis_tbl,
         clust_clinic$variables,
         ~map_lgl(.x[.y], is.numeric)) %>% 
    map(~ifelse(.x, 'years', '% of cluster'))
  
# Testing ----
  
  insert_msg('Testing')
  
  ## test results
  
  clust_clinic$test_results <- 
    list(data = clust_clinic$analysis_tbl, 
         variables = clust_clinic$variables, 
         types = clust_clinic$types) %>% 
    future_pmap(function(data, variables, types) compare_variables(data, 
                                                                   variables = variables, 
                                                                   split_factor = 'clust_id', 
                                                                   what = 'eff_size', 
                                                                   types = types, 
                                                                   exact = FALSE, 
                                                                   ci = FALSE, 
                                                                   pub_styled = TRUE, 
                                                                   adj_method = 'BH'), 
                .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## significant and near-significant (p < 0.01) factors
  
  clust_clinic$top_factors <- clust_clinic$test_results %>% 
    map(filter, p_adjusted < 0.1)
  
# Plotting ------
  
  insert_msg('Plotting')
  
  clust_clinic$plots <- 
    list(data = clust_clinic$analysis_tbl, 
         cohort = globals$cohort_lexicon$label, 
         variable = clust_clinic$variables, 
         plot_type = clust_clinic$plot_types, 
         stats = clust_clinic$test_results, 
         plot_title = map(clust_clinic$variables, 
                          exchange, 
                          dict = globals$var_lexicon, 
                          key = 'variable', 
                          value = 'label'), 
         y_labs = clust_clinic$y_labels) %>% 
    pmap(function(data, 
                  cohort, 
                  variable, 
                  plot_type, 
                  stats, 
                  plot_title, 
                  y_labs) list(variable = variable, 
                               type = plot_type, 
                               plot_title = paste(plot_title, 
                                                  cohort, 
                                                  sep = ', '), 
                               plot_subtitle = stats$plot_cap, 
                               y_lab = y_labs) %>% 
           pmap(plot_variable, 
                data, 
                split_factor = 'clust_id', 
                cust_theme = globals$common_theme, 
                x_lab = 'Lymphokine cluster', 
                scale = 'percent') %>% 
           map(~.x + 
                 scale_fill_brewer() + 
                 scale_x_discrete(labels = .x$labels$tag %>% 
                                    stri_split(fixed = '\n') %>% 
                                    map(stri_replace, 
                                        fixed = ': ', 
                                        replacement = '\n')) + 
                 theme(plot.tag = element_blank())) %>% 
           set_names(variable))
  
# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()