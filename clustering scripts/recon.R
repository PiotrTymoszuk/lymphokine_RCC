# Differences in signatures of the Recon2 metabolism subsystems between the 
# lymphokine clusters

  insert_head()
  
# container ------
  
  clust_recon <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  clust_recon$variables <- sig$recon_lexicon$variable %>% 
    unique
  
  clust_recon$variables <- 
    clust_recon$variables[clust_recon$variables != 'unassigned']
  
  clust_recon$var_lexicon <- sig$recon_lexicon %>% 
    filter(!duplicated(variable))
  
  ## analysis tables
  
  clust_recon$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('clust_recon$analysis_tbl$', i, 
                           ' <- ', i, '$recon')))
    
  }
  
  clust_recon$analysis_tbl <- gene_clust$clust_obj %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(., clust_recon$analysis_tbl, 
         left_join, by = 'patient_id')
  
  clust_recon$variables <- clust_recon$analysis_tbl %>% 
    map(~clust_recon$variables[clust_recon$variables %in% names(.x)])
  
# Serial testing -----
  
  insert_msg('Serial testing')
  
  ## two-tailed T test
  
  clust_recon$test_results <- 
    list(data = clust_recon$analysis_tbl, 
         variables = clust_recon$variables ) %>% 
    future_pmap(test_two_groups, 
                split_fct = 'clust_id', 
                adj_method = 'BH', 
                .parallel = FALSE) %>% 
    map(mutate, 
        sign_label = exchange(response, 
                              dict = clust_recon$var_lexicon, 
                              key = 'variable', 
                              value = 'label'))
  
  ## significantly regulated signatures
  
  clust_recon$significant <- clust_recon$test_results %>% 
    map(filter, significant == 'yes')
  
# Volcano plots with the results ------
  
  insert_msg('Volcano plots')
  
  clust_recon$volcano_plots <- 
    list(data = clust_recon$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         label_variable = 'sign_label', 
         top_significant = 20, 
         x_lab = expression(Delta * ' ssGSEA score, chemox high vs chemox low'), 
         y_lab = expression('-log'[10] * ' pFDR'), 
         txt_size = 2.4, 
         label_type = 'text', 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          geom_vline(xintercept = 0, 
                     linetype = 'dashed') + 
          labs(subtitle = .x$labels$tag) + 
          theme(plot.tag = element_blank()))
  
# Top regulated plots ------
  
  insert_msg('Top regulated plots')
  
  clust_recon$top_plots <- 
    list(data = clust_recon$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_top, 
         regulation_variable = 'estimate', 
         label_variable = 'sign_label', 
         p_variable = 'p_adjusted', 
         lower_ci_variable = 'lower_ci', 
         upper_ci_variable = 'upper_ci', 
         top_regulated = 20, 
         plot_subtitle = 'Top regulated Recon2 metabolic subsystem signatures', 
         x_lab = expression(Delta * ' ssGSEA score, chemox high vs chemox low'), 
         cust_theme = globals$common_theme)
  
# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()