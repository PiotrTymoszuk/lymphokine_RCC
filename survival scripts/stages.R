# Analysis of overall survival in low-stage (T1 - T2) and  high-stage cancers
# Uni-variable Cox modeling
#
# Done for the TCGA, e-MTAB-1980 and RECA-EU collectives with first-order and 
# expression spline models

  insert_head()
  
# container ------
  
  stage_cox <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## analysis cohorts
  
  stage_cox$cohorts <- c('tcga', 'emtab1980', 'reca')
  
  stage_cox$cohort_labs <- stage_cox$cohorts %>% 
    exchange(dict = globals$cohort_lexicon, 
             key = 'cohort',
             value = 'label') %>% 
    map(~paste(.x, c('T1 - T2', 'T3 - T4'))) %>% 
    unlist
  
  ## analysis tables with normalized expression
  ## of the genes of interest
  
  stage_cox$analysis_tbl <- list()
  
  for(i in stage_cox$cohorts) {
    
    eval(parse_expr(paste0('stage_cox$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  stage_cox$analysis_tbl <- stage_cox$analysis_tbl %>% 
    map(select, 
        patient_id, 
        pt_stage, 
        os_days,
        death, 
        any_of('tissue_type'), 
        all_of(globals$cxcl_genes))
  
  stage_cox$analysis_tbl$tcga <- stage_cox$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  for(i in globals$cxcl_genes) {
    
    stage_cox$analysis_tbl <- stage_cox$analysis_tbl %>% 
      map(mutate, 
          !!i := scale(.data[[i]])[, 1])
    
  }
  
  ## splitting by the stage
  
  stage_cox$analysis_tbl <- stage_cox$analysis_tbl %>% 
    map(mutate, 
        pt_stage = stri_extract(pt_stage, regex = 'T\\d{1}'), 
        progression = ifelse(pt_stage %in% c('T1', 'T2'), 
                             'early', 'late'), 
        progression = factor(progression, c('early', 'late'))) %>% 
    map(~filter(.x, complete.cases(.x)))
  
  stage_cox$analysis_tbl <- stage_cox$analysis_tbl %>% 
    map(dlply, 'progression', as_tibble) %>% 
    unlist(recursive = FALSE)
  
# Construction of the models ------
  
  insert_msg('Model construction')
  
  ## no spline formulas
  
  stage_cox$formulas <- globals$cxcl_genes %>% 
    set_names(globals$cxcl_genes) %>% 
    map(~paste('Surv(os_days, death) ~', .x)) %>% 
    map(as.formula)
  
  ## no spline models, OS
  
  for(i in names(stage_cox$formulas)) {
    
    stage_cox$models[[i]] <- stage_cox$analysis_tbl %>% 
      map(~call2(.fn = 'coxph', 
                 formula = stage_cox$formulas[[i]], 
                 data = .x)) %>% 
      map(eval) %>% 
      map2(., stage_cox$analysis_tbl, as_coxex)
    
  }
  
  ## spline formulas
  
  stage_cox$spline_formulas <- globals$cxcl_genes %>% 
    set_names(globals$cxcl_genes) %>% 
    map(~paste0('Surv(os_days, death) ~ pspline(', .x, ', df = 4)')) %>% 
    map(as.formula)
  
  ## spline models, OS
  
  for(i in names(stage_cox$spline_formulas)) {
    
    stage_cox$spline_models[[i]] <- stage_cox$analysis_tbl %>% 
      map(~call2(.fn = 'coxph', 
                 formula = stage_cox$spline_formulas[[i]], 
                 data = .x)) %>% 
      map(eval) %>% 
      map2(., stage_cox$analysis_tbl, as_coxex)
    
  }
  
  
  
# Model assumptions --------
  
  insert_msg('Model assumptions')
  
  ## proportional hazard assumption
  
  stage_cox$assumptions <- stage_cox$models %>% 
    map(~map(.x, summary, 'assumptions')) %>% 
    map(compress, names_to = 'cohort')
  
  stage_cox$spline_assumptions <- stage_cox$spline_models %>% 
    map(~map(.x, summary, 'assumptions')) %>% 
    map(compress, names_to = 'cohort')
  
# Fit stats --------
  
  insert_msg('Fit stats')
  
  ## no spline models
  
  stage_cox$fit <- stage_cox$models %>% 
    map(~map(.x, summary, 'fit')) %>% 
    map(compress, names_to = 'cohort') %>% 
    map(mutate, 
        cohort_lab = stage_cox$cohort_labs, 
        n_lab = paste0('total: n = ', n_complete, 
                       ', events: n = ', n_events))
  
  ## spline models
  
  stage_cox$spline_fit <- stage_cox$spline_models %>% 
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
  
  stage_cox$inference <- stage_cox$models %>% 
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
  
  stage_cox$forest_titles <- globals$cxcl_genes %>% 
    map(~paste('<b><em>', .x, '</em>, overall survival</b>'))
  
  stage_cox$forest_cohorts <- stage_cox$fit %>% 
    map(mutate, 
        cohort_lab = paste(cohort_lab, n_lab, sep = '\n'))
  
  stage_cox$forest_cohorts <- stage_cox$forest_cohorts %>% 
    map(~set_names(.x$cohort_lab, .x$cohort))
  
  ## Forest plots
  
  stage_cox$forest <- 
    list(data = stage_cox$inference, 
         plot_title = stage_cox$forest_titles) %>% 
    pmap(plot_top, 
         label_variable = 'cohort', 
         p_variable = 'p_value', 
         regulation_variable = 'estimate', 
         lower_ci_variable = 'lower_ci', 
         upper_ci_variable = 'upper_ci', 
         plot_subtitle = 'Uni-variable Cox modeling, linear term', 
         x_lab = 'log HR, 95% CI', 
         cust_theme = globals$common_theme) %>% 
    map2(., stage_cox$forest_cohorts, 
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
  
  stage_cox$spline_inference <- stage_cox$spline_models %>% 
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
  
  stage_cox$spline_titles <- globals$cxcl_genes %>% 
    map(function(gene) stage_cox$cohort_labs %>% 
          map(~paste0('<b><em>', gene, '</em>, ', .x, '</b>'))) %>% 
    set_names(globals$cxcl_genes)
  
  ## plotting
  
  stage_cox$spline_plots <- list(model = stage_cox$spline_models, 
                               title = stage_cox$spline_titles, 
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
  
# END ------
  
  rm(i)
  
  insert_tail()