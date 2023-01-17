# Uni-variable Cox modeling of overall survival
#
# Survival data are available for TCGA, RECA-EU, E-MTAB-1980 and CheckMate 
# cohorts. Models with and without spline of expression of the genes of interest
# are constructed

  insert_head()
  
# container -------
  
  uni_cox <- list()
  
# Analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## analysis tables, normalized expression values
  
  uni_cox$analysis_tbl <- list()
  
  for(i in globals$os_cohorts) {
    
    eval(parse_expr(paste0('uni_cox$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  uni_cox$analysis_tbl <- uni_cox$analysis_tbl %>% 
    map(select, 
        patient_id, 
        os_days,
        death, 
        any_of('tissue_type'), 
        all_of(globals$cxcl_genes)) %>% 
    map(~filter(.x, complete.cases(.x)))
  
  uni_cox$analysis_tbl$tcga <- 
    uni_cox$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  for(i in globals$cxcl_genes) {
    
    uni_cox$analysis_tbl <- uni_cox$analysis_tbl %>% 
      map(mutate, 
          !!i := scale(.data[[i]])[, 1])
    
  }
  
# Construction of the models ------
  
  insert_msg('Model construction')
  
  ## no spline formulas
  
  uni_cox$formulas <- globals$cxcl_genes %>% 
    set_names(globals$cxcl_genes) %>% 
    map(~paste('Surv(os_days, death) ~', .x)) %>% 
    map(as.formula)
  
  ## no spline models, OS
  
  for(i in names(uni_cox$formulas)) {
    
    uni_cox$models[[i]] <- uni_cox$analysis_tbl %>% 
      map(~call2(.fn = 'coxph', 
                 formula = uni_cox$formulas[[i]], 
                 data = .x)) %>% 
      map(eval) %>% 
      map2(., uni_cox$analysis_tbl, as_coxex)
    
  }
  
  ## spline formulas
  
  uni_cox$spline_formulas <- globals$cxcl_genes %>% 
    set_names(globals$cxcl_genes) %>% 
    map(~paste0('Surv(os_days, death) ~ pspline(', .x, ', df = 4)')) %>% 
    map(as.formula)

  ## spline models, OS
  
  for(i in names(uni_cox$spline_formulas)) {
    
    uni_cox$spline_models[[i]] <- uni_cox$analysis_tbl %>% 
      map(~call2(.fn = 'coxph', 
                 formula = uni_cox$spline_formulas[[i]], 
                 data = .x)) %>% 
      map(eval) %>% 
      map2(., uni_cox$analysis_tbl, as_coxex)
    
  }
  
# Model assumptions --------
  
  insert_msg('Model assumptions')
  
  ## proportional hazard assumption
  
  uni_cox$assumptions <- uni_cox$models %>% 
    map(~map(.x, summary, 'assumptions')) %>% 
    map(compress, names_to = 'cohort')
    
  uni_cox$spline_assumptions <- uni_cox$spline_models %>% 
    map(~map(.x, summary, 'assumptions')) %>% 
    map(compress, names_to = 'cohort')
  
# Fit stats --------
  
  insert_msg('Fit stats')
  
  ## no spline models
  
  uni_cox$fit <- uni_cox$models %>% 
    map(~map(.x, summary, 'fit')) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, 
        cohort_lab = exchange(cohort, 
                              dict = globals$cohort_lexicon, 
                              key = 'cohort', 
                              value = 'label'), 
        n_lab = paste0('total: n = ', n_complete, 
                       ', events: n = ', n_events))
  
  ## spline models
  
  uni_cox$spline_fit <- uni_cox$spline_models %>% 
    map(~map(.x, as_coxph) %>% 
          map(~tibble(aic = AIC(.x), 
                      bic = BIC(.x), 
                      raw_rsq = survMisc::rsq(.x)[['mev']], 
                      c_index = concordance(.x)[['concordance']], 
                      c_index_se = concordance(.x)[['var']]^0.5)) %>% 
          map(mutate, 
              lower_ci = c_index + qnorm(0.025) * c_index_se, 
              upper_ci = c_index + qnorm(0.975) * c_index_se)) %>% 
    map(compress, names_to = 'cohort')
  
# Inference for the no-spline models -------
  
  insert_msg('Inference, no spline models and Forest plots')
  
  uni_cox$inference <- uni_cox$models %>% 
    map(~map(.x, summary, 'inference')) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, 
        or = exp(estimate), 
        or_lower = exp(lower_ci), 
        or_upper = exp(upper_ci), 
        plot_cap = paste0(signif(or, 2), 
                          ' [', signif(or_lower, 2), ' - ', 
                          signif(or_upper, 2), ']'))
  
  ## Forest plot titles and x labels
  
  uni_cox$forest_titles <- globals$cxcl_genes %>% 
    map(~paste0('<b><em>', .x, '</em>, overall survival</b>'))
  
  uni_cox$forest_cohorts <- uni_cox$fit %>% 
    map(mutate, 
        cohort_lab = paste(cohort_lab, n_lab, sep = '\n'))
  
  uni_cox$forest_cohorts <- uni_cox$forest_cohorts %>% 
    map(~set_names(.x$cohort_lab, .x$cohort))
  
  ## Forest plots
  
  uni_cox$forest <- 
    list(data = uni_cox$inference, 
         plot_title = uni_cox$forest_titles) %>% 
    pmap(plot_top, 
         label_variable = 'cohort', 
         p_variable = 'p_value', 
         regulation_variable = 'estimate', 
         lower_ci_variable = 'lower_ci', 
         upper_ci_variable = 'upper_ci', 
         plot_subtitle = 'Uni-variable Cox modeling, linear term', 
         x_lab = 'log HR, 95% CI', 
         cust_theme = globals$common_theme) %>% 
    map2(., uni_cox$forest_cohorts, 
         ~.x + 
          geom_text(aes(label = plot_cap), 
                    size = 2.6, 
                    hjust = 0.5, 
                    vjust = -0.8) + 
          scale_y_discrete(labels = .y) + 
          theme(plot.title = element_markdown()))
  
# Inference for the spline models ----
  
  insert_msg('Inference for the spline models')
  
  ## inference: data for nodes can be obtained
  ## by calling `summary()`. I'm extracting the survival (log OR)
  ## predictions for the terms. See:
  ## https://cran.r-project.org/web/packages/survival/vignettes/splines.pdf
  
  uni_cox$spline_inference <- uni_cox$spline_models %>% 
    map(~map(.x, as_coxph) %>% 
          map(termplot, 
              plot = FALSE, 
              se = TRUE) %>% 
          map(~.x[[1]]) %>% 
          map(transmute, 
              expression = x, 
              estimate = y, 
              se = se, 
              lower_ci = estimate + qnorm(0.025) * se, 
              upper_ci = estimate + qnorm(0.975) * se) %>% 
          map(as_tibble))
  
  ## plot titles
  
  uni_cox$spline_titles <- globals$cxcl_genes %>% 
    map(function(gene) globals$cohort_lexicon %>% 
          filter(cohort %in% globals$os_cohorts) %>% 
          .$label %>% 
          map(~paste0('<b><em>', gene, '</em>, ', .x, '</b>'))) %>% 
    set_names(globals$cxcl_genes)
  
  ## plotting
  
  uni_cox$spline_plots <- list(model = uni_cox$spline_models, 
                               title = uni_cox$spline_titles, 
                               color = c('coral4', 
                                         'coral3', 
                                         'coral2', 
                                         'steelblue3', 
                                         'darkolivegreen4', 
                                         'gray60', 
                                         'firebrick3', 
                                         'darkorange3')) %>% 
    pmap(function(model, title, color) list(coxex_model = model, 
                                            plot_title = title) %>% 
           pmap(plot_spline_cox, 
                color = color) %>% 
           map(~.x +
                 geom_rug(sides = 'b', 
                          color = 'gray60') + 
                 theme(plot.title = element_markdown())))

# significance of the non-linear terms -----
  
  insert_msg('Significance of the on-linear gene expression terms')
  
  uni_cox$term_significance <- uni_cox$spline_models %>% 
    map(~map(.x, as_coxph) %>% 
          map(summary) %>%
          map(~.x$coefficients) %>% 
          map(as.data.frame) %>% 
          map(mutate, term = c('linear', 'spline')) %>% 
          map(as_tibble) %>% 
          compress(names_to = 'cohort'))
  
# Plots of the p values for the spline terms ------
  
  insert_msg('Plots of the p values of the spline terms')
  
  uni_cox$spline_significance_plots <- 
    list(data = map(uni_cox$term_significance, 
                    filter, term == 'spline'), 
       plot_title = uni_cox$forest_titles) %>% 
    pmap(plot_signifcant, 
         p_variable = 'p', 
         label_variable = 'cohort', 
         plot_subtitle = 'Significance of spline terms', 
         cust_theme = globals$common_theme) %>% 
    map2(., uni_cox$forest_cohorts, 
         ~.x + 
           scale_y_discrete(labels = .y) + 
           theme(plot.title = element_markdown()))
  
# END ------
  
  rm(i)
  
  insert_tail()