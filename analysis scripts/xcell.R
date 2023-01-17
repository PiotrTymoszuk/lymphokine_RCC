# Correlation of the CXCL9/10/11 and CXCL11 with immune cell infiltration.
# Spearman correlation.

  insert_head()

# container list ----

  xcell <- list()

# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# globals -----

  insert_msg('Globals setup')

  ## analysis table

  xcell$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('xcell$analysis_tbl$xcell$', i, 
                           ' <- ', i, '$xcell')))
    
    eval(parse_expr(paste0('xcell$analysis_tbl$expression$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  xcell$analysis_tbl$expression[c('tcga', 'gse167093')] <- 
    xcell$analysis_tbl$expression[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  xcell$analysis_tbl$expression <- 
    xcell$analysis_tbl$expression %>% 
    map(~.x[c('patient_id', globals$cxcl_genes)])
  
  xcell$analysis_tbl <- 
    map2(xcell$analysis_tbl$xcell, 
         xcell$analysis_tbl$expression, 
         inner_join, by = 'patient_id')
    
  ## variable pairs

  xcell$pairs <- globals$cxcl_genes %>%
    map(function(cxc) map(globals$xcell_lexicon$variable,
                          ~c(cxc, .x)) %>%
          set_names(c(globals$xcell_lexicon$variable))) %>%
    set_names(globals$cxcl_genes)

# correlations ----

  insert_msg('Serial correlation')

  xcell$test_results <- set_names(names(xcell$analysis_tbl),
                                  names(xcell$analysis_tbl)) %>%
    future_map(function(cohort) xcell$pairs %>%
                 map(function(cxc) map_dfr(cxc,
                                           ~correlate_variables(xcell$analysis_tbl[[cohort]],
                                                                variables = .x,
                                                                what = 'correlation',
                                                                type = 'spearman',
                                                                pub_styled = FALSE,
                                                                simplify_p = FALSE))) %>%
                 map_dfr(mutate, 
                         p_adjusted = p.adjust(p_value, 'BH')), 
               .options = furrr_options(seed = TRUE))
  
  xcell$test_results <- xcell$test_results %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 
                             'bold', 'plain'),
        plot_lab = signif(estimate, 2), 
        estimate_recoded = ifelse(p_adjusted < 0.05, 
                                  estimate, NA))

# Summary bubble plots -----

  insert_msg('Summary plots')

  xcell$bubble <- list(correlation_data = xcell$test_results, 
                       plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_immune_bubble, 
         variable_lexicon = globals$xcell_lexicon)
  
# Summary simplified bubble plots --------
  
  insert_msg('Summary simplified bubble plots')
  
  ## no labels, non-significant correlations are squished to 0
  
  xcell$heat_bubble <- 
    list(correlation_data = xcell$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_simplified_bubble, 
         variable_lexicon = globals$xcell_lexicon)
  
# Point plots for the CD8 and CXCL9/CXCR3 correlations -------

  insert_msg('Correlation point plots')

  ## plot captions

  xcell$point_captions <- xcell$test_results %>% 
    map(mutate, 
        plot_cap = paste0(signif(estimate, 2), 
                          ' [', signif(lower_ci, 2), 
                          ' - ', signif(upper_ci, 2), ']'), 
        plot_cap = paste(plot_cap, significance, sep = ', '), 
        plot_cap = paste(plot_cap, n, sep = ', n = ')) %>% 
    map(select, 
        variable1, 
        variable2, 
        plot_cap)
  
  ## plot titles and x, y axis labels
  
  xcell$point_titles <- 
    map2(xcell$point_captions, 
         globals$cohort_lexicon$label, 
         function(var, cohort) var$variable1 %>% 
           map(~paste0('<b><em>', .x, '</em>, ', cohort, '</b>')))
  
  xcell$point_x_labs <- xcell$point_captions %>% 
    map(function(var) var$variable1 %>% 
           map(~paste0('<em>', .x, '</em>, log<sub>2</sub> expression')))
  
  xcell$point_y_labs <- xcell$point_captions %>% 
    map(function(var) var$variable2 %>% 
          map(~paste0(.x, ', xCell')))

  ## point plots
  
  xcell$point_plots <- list(data = xcell$analysis_tbl, 
                            stats = xcell$point_captions, 
                            title = xcell$point_titles, 
                            x_lab = xcell$point_x_labs,
                            y_lab = xcell$point_y_labs) %>% 
    future_pmap(function(data, stats, title, x_lab, y_lab) list(x = stats$variable1, 
                                                                y = stats$variable2, 
                                                                z = stats$plot_cap, 
                                                                v = title, 
                                                                w = x_lab, 
                                                                q = y_lab) %>% 
                  pmap(function(x, y, z, v, w, q) plot_correlation(data = data, 
                                                                   variables = c(x, y), 
                                                                   type = 'correlation', 
                                                                   point_hjitter = 0, 
                                                                   point_wjitter = 0, 
                                                                   cust_theme = globals$common_theme, 
                                                                   plot_title = v, 
                                                                   plot_subtitle = z, 
                                                                   x_lab = w, 
                                                                   y_lab = q, 
                                                                   show_trend = FALSE) + 
                         geom_smooth(method = 'gam') + 
                         scale_y_continuous(trans = 'pseudo_log') + 
                         theme(plot.title = element_markdown(), 
                               axis.title.x = element_markdown(), 
                               plot.tag = element_blank())) %>% 
                  set_names(paste(stats$variable1, 
                                  stats$variable2, 
                                  sep = ':')), 
                .options = furrr_options(seed = TRUE))
  
# END ------
  
  rm(i)
  
  plan('sequential')

  insert_tail()
