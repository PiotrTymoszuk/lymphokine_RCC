# Correlation of the CXCL9/10/11 and CXCL11 with immune cell infiltration.
# Spearman correlation.

  insert_head()

# container list ----

  quantiseq <- list()

# globals -----

  insert_msg('Globals setup')

  ## analysis table

  quantiseq$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('quantiseq$analysis_tbl$quantiseq$', i, 
                           ' <- ', i, '$quantiseq')))
    
    eval(parse_expr(paste0('quantiseq$analysis_tbl$expression$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  quantiseq$analysis_tbl$expression[c('tcga', 'gse167093')] <- 
    quantiseq$analysis_tbl$expression[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  quantiseq$analysis_tbl$expression <- 
    quantiseq$analysis_tbl$expression %>% 
    map(~.x[c('patient_id', globals$cxcl_genes)])
  
  quantiseq$analysis_tbl <- 
    map2(quantiseq$analysis_tbl$quantiseq, 
         quantiseq$analysis_tbl$expression, 
         inner_join, by = 'patient_id')
    
  ## variable pairs

  quantiseq$pairs <- globals$cxcl_genes %>%
    map(function(cxc) map(globals$quantiseq_lexicon$variable,
                          ~c(cxc, .x)) %>%
          set_names(c(globals$quantiseq_lexicon$variable))) %>%
    set_names(globals$cxcl_genes)

# correlations ----

  insert_msg('Serial correlation')

  quantiseq$test_results <- set_names(names(quantiseq$analysis_tbl),
                                       names(quantiseq$analysis_tbl)) %>%
    map(function(cohort) quantiseq$pairs %>%
          map(function(cxc) map_dfr(cxc,
                                    ~correlate_variables(quantiseq$analysis_tbl[[cohort]],
                                                         variables = .x,
                                                         what = 'correlation',
                                                         type = 'spearman',
                                                         pub_styled = FALSE,
                                                         simplify_p = FALSE))) %>%
          map_dfr(mutate, p_adjusted = p.adjust(p_value, 'BH')))

  quantiseq$test_results <- quantiseq$test_results %>%
    map(mutate,
        significant = ifelse(p_adjusted < 0.05, 
                             'bold', 'plain'),
        plot_lab = signif(estimate, 2), 
        estimate_recoded = ifelse(p_adjusted < 0.05, 
                                  estimate, NA))

# Summary bubble plots -----

  insert_msg('Summary plots')

  quantiseq$bubble <- list(correlation_data = quantiseq$test_results, 
                           plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_immune_bubble, 
         variable_lexicon = globals$quantiseq_lexicon)
  
# Summary simplified bubble plots --------
  
  insert_msg('Summary simplified bubble plots')
  
  ## no labels, non-significant correlations are squished to 0
  
  quantiseq$heat_bubble <- 
    list(correlation_data = quantiseq$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_simplified_bubble, 
         variable_lexicon = globals$quantiseq_lexicon, 
         max_size = 5)

# Point plots for the CD8 and CXCL9/CXCR3 correlations -------

  insert_msg('Correlation point plots')

  ## plot captions

  quantiseq$point_captions <- quantiseq$test_results %>% 
    map(filter, variable2 == 'T cell CD8+') %>% 
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
  
  ## plot titles and x axis labels
  
  quantiseq$point_titles <- 
    map2(quantiseq$point_captions, 
         globals$cohort_lexicon$label, 
         function(var, cohort) var$variable1 %>% 
           map(~paste0('<b><em>', .x, '</em>, ', cohort, '</b>')))
  
  quantiseq$point_x_labs <- quantiseq$point_captions %>% 
    map(function(var) var$variable1 %>% 
           map(~paste0('<em>', .x, '</em>, log<sub>2</sub> expression')))

  ## point plots
  
  quantiseq$point_plots <- list(data = quantiseq$analysis_tbl, 
                                 stats = quantiseq$point_captions, 
                                 title = quantiseq$point_titles, 
                                 x_lab = quantiseq$point_x_labs) %>% 
    pmap(function(data, stats, title, x_lab) list(x = stats$variable1, 
                                                  y = stats$variable2, 
                                                  z = stats$plot_cap, 
                                                  v = title, 
                                                  w = x_lab) %>% 
           pmap(function(x, y, z, v, w) plot_correlation(data = data, 
                                                         variables = c(x, y), 
                                                         type = 'correlation', 
                                                         point_hjitter = 0, 
                                                         point_wjitter = 0, 
                                                         cust_theme = globals$common_theme, 
                                                         plot_title = v, 
                                                         plot_subtitle = z, 
                                                         x_lab = w, 
                                                         y_lab = 'CD8+ T cells, QuanTIseq', 
                                                         show_trend = FALSE) + 
                  geom_smooth(method = 'gam') + 
                  scale_y_continuous(trans = 'pseudo_log') + 
                  theme(plot.title = element_markdown(), 
                        axis.title.x = element_markdown(), 
                        plot.tag = element_blank())) %>% 
           set_names(stats$variable1))
  
# END ------
  
  rm(i)

  insert_tail()
