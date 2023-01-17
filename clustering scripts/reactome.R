# Differences in signatures of the Recon2 metabolism subsystems between the 
# lymphokine clusters

  insert_head()
  
# container ------
  
  clust_reactome <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  clust_reactome$variables <- sig$reactome_lexicon$variable %>% 
    unique

  clust_reactome$var_lexicon <- sig$reactome_lexicon %>% 
    filter(!duplicated(variable))
  
  ## analysis tables
  
  clust_reactome$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('clust_reactome$analysis_tbl$', i, 
                           ' <- ', i, '$reactome')))
    
  }
  
  clust_reactome$analysis_tbl <- gene_clust$clust_obj %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(., clust_reactome$analysis_tbl, 
         left_join, by = 'patient_id')
  
  clust_reactome$variables <- clust_reactome$analysis_tbl %>% 
    map(~clust_reactome$variables[clust_reactome$variables %in% names(.x)])
  
# Serial testing -----
  
  insert_msg('Serial testing')
  
  ## two-tailed T test
  
  clust_reactome$test_results <- 
    list(data = clust_reactome$analysis_tbl, 
         variables = clust_reactome$variables ) %>% 
    future_pmap(test_two_groups, 
                split_fct = 'clust_id', 
                adj_method = 'BH', 
                .parallel = FALSE) %>% 
    map(mutate, 
        sign_label = exchange(response, 
                              dict = clust_reactome$var_lexicon, 
                              key = 'variable', 
                              value = 'label'))
  
  ## significantly regulated signatures
  
  clust_reactome$significant <- clust_reactome$test_results %>% 
    map(filter, significant == 'yes')
  
# Identification of the common up- and downregulated signatures -----
  
  insert_msg('Common signatures')
  
  ## general, significant signatures shared by all datasets
  
  clust_reactome$common_signatures <- clust_reactome$significant %>% 
    map(~.x$response) %>% 
    reduce(intersect)
  
  ## signatures shared by at least 6 out of 8 datasets
  
  clust_reactome$cmm_sets <- utils::combn(names(clust_reactome$significant),
                                          m = 7,
                                          simplify = FALSE)
  
  clust_reactome$shared_signatures <- clust_reactome$cmm_sets %>% 
    map(~clust_reactome$significant[.x]) %>% 
    map(~map(.x, ~.x$response)) %>% 
    map(reduce, intersect) %>% 
    reduce(union)
  
# Volcano plots with the results ------
  
  insert_msg('Volcano plots')
  
  clust_reactome$volcano_plots <- 
    list(data = clust_reactome$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         label_variable = 'sign_label', 
         top_significant = 10, 
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
  
  clust_reactome$top_plots <- 
    list(data = clust_reactome$significant, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_top, 
         regulation_variable = 'estimate', 
         label_variable = 'sign_label', 
         p_variable = 'p_adjusted', 
         lower_ci_variable = 'lower_ci', 
         upper_ci_variable = 'upper_ci', 
         top_regulated = 20, 
         plot_subtitle = 'Top regulated Reactome signatures', 
         x_lab = expression(Delta * ' ssGSEA score, chemox high vs chemox low'), 
         cust_theme = globals$common_theme)
  
# Bubble plot with the significantly regulated common signatures -----
  
  insert_msg('Bubble plot')
  
  ## signatures shared by all datasets are plotted
  
  clust_reactome$bubble_plot <- clust_reactome$significant %>% 
    map(filter, response %in% clust_reactome$common_signatures) %>% 
    map(select, 
        response, sign_label, estimate) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(regulation = ifelse(estimate > 0, 
                               'upregulated', 'downregulated'), 
           cohort = factor(cohort, globals$cohort_lexicon$cohort)) %>% 
    ggplot(aes(x = cohort, 
               y = reorder(sign_label, estimate), 
               fill = estimate, 
               size = estimate)) + 
    geom_point(shape = 21) + 
    scale_size_area(name = '\u0394 ssGSEA score\nchemox high vs chemox low') + 
    scale_x_discrete(labels = globals$cohort_lexicon$label) + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'white', 
                         high = 'firebrick', 
                         midpoint = 0,
                         name = '\u0394 ssGSEA score\nchemox high vs chemox low') + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
    labs(title = 'Commmon significantly regulated Reactome pathways')

# Volcano plots with the common regulated processes labeled ------
  
  insert_msg('Volcano plots, common processes labeled')
  
  ## plotting tables
  
  clust_reactome$cmm_volcano_tbl <- clust_reactome$test_results %>% 
    map(filter, response %in% clust_reactome$common_signatures) %>% 
    map(mutate, 
        sign_label = reactome_labeller(sign_label), 
        color = ifelse(sign_label %in% c('PD1 signaling', 
                                         'IFN gamma signaling', 
                                         'IFN alpha beta signaling', 
                                         'IFN signaling'), 
                      'bold', 'plain'), 
        regulation = ifelse(estimate > 0, 'upregulated', 'downregulated'))
  
  ## volcano plots
  
  clust_reactome$cmm_volcano_plots <- 
    list(data = clust_reactome$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         x_lab = expression(Delta * ' ssGSEA score, chemox high vs chemox low'), 
         y_lab = expression('-log'[10] * ' pFDR'), 
         txt_size = 2.4, 
         point_alpha = 0.15, 
         label_type = 'text', 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          geom_vline(xintercept = 0, 
                     linetype = 'dashed') + 
          labs(subtitle = .x$labels$tag) + 
          theme(plot.tag = element_blank()))
    
  ## overlay with the common regulated signatures
  
  clust_reactome$cmm_volcano_plots <- 
    map2(clust_reactome$cmm_volcano_plots, 
         clust_reactome$cmm_volcano_tbl, 
         ~.x + 
           geom_point(data = .y, 
                      shape = 21, 
                      size = 2) + 
           geom_text_repel(data = .y, 
                           aes(label = sign_label, 
                               color = color), 
                           size = 2.5) + 
           scale_color_manual(values = c(bold = 'black', 
                                         plain = 'gray40')) + 
           guides(color = 'none'))
  
# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()