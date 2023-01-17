# Multi-parameter overall survival modeling.
#
# Explanatory variables: lymphokine genes of interest
# modeling with Ridge Cox regression (maximal regularization of the model)
#
# The TCGA cohort serves as a training collective. The linear predictor scores
# of the TCGA GLMNET model are predicted for the remaining collectives 
# and the association of the linear predictor scores with survival is modeled
# by uni-variable Cox regression

  insert_head()
  
# container ------
  
  multi_cox <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## genes of interest
  
  multi_cox$genes <- globals$cxcl_genes
  
  ## analysis tables
  ## adding a day to OS, because the algorithm does not tolerate
  ## immediate censoring - no effects on modeling results expected
  
  multi_cox$analysis_tbl <- list()
  
  for(i in globals$os_cohorts) {
    
    eval(parse_expr(paste0('multi_cox$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  multi_cox$analysis_tbl <- multi_cox$analysis_tbl %>% 
    map(select, 
        patient_id, 
        os_days,
        death, 
        any_of('tissue_type'), 
        all_of(multi_cox$genes)) %>% 
    map(~filter(.x, complete.cases(.x)))
  
  multi_cox$analysis_tbl$tcga <- 
    multi_cox$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  multi_cox$analysis_tbl <- multi_cox$analysis_tbl %>% 
    map(mutate, 
        os_days = os_days + 1)
  
  for(i in multi_cox$genes) {
    
    multi_cox$analysis_tbl <- multi_cox$analysis_tbl %>% 
      map(mutate, 
          !!i := scale(.data[[i]])[, 1])
    
  }
  
  ## survival objects and explanatory variable matrices
  ## no higher order terms included in the models
  ## since no significant effects for splines were detected
  ## in uni-variable Cox modeling
  
  multi_cox$y <- multi_cox$analysis_tbl %>% 
    map(~Surv(.x$os_days, .x$death))
  
  multi_cox$x <- multi_cox$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(~.x[multi_cox$genes]) %>% 
    map(~model.matrix(~., data = .x))
  
  ## alpha and CV folds
  
  multi_cox$alpha <- 0 
  multi_cox$n_rep <- 200
  
  ## folds for CV: TCGA is the training cohort
  
  set.seed(1234)
  
  multi_cox$folds <- 1:multi_cox$n_rep %>% 
    map(~createFolds(y = multi_cox$analysis_tbl$tcga$patient_id, 
                     k = 10, 
                     list = FALSE, 
                     returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:multi_cox$n_rep))
  
# parallel backend -----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# lambda tuning -------
  
  insert_msg('Lambda tuning, TCGA')
  
  ## model construction
  
  multi_cox$lambda_tune <- multi_cox$folds %>% 
    future_map(~cv.glmnet(x = multi_cox$x$tcga, 
                          y = multi_cox$y$tcga, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = multi_cox$alpha), 
               .options = furrr_options(seed = TRUE))
  
  ## extracting the optimal lambdas
  
  multi_cox$lambda_tbl <- multi_cox$lambda_tune %>% 
    map(~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo', 
                        'lambda.min', 'lambda.1se')])) %>% 
    map_dfr(filter, lambda == lambda.1se)
  
  multi_cox$opt_lambda <- multi_cox$lambda_tbl %>% 
    filter(cvm == min(cvm))
  
# Building the GLMNET model -------
  
  insert_msg('Building the GLMNET model')

  multi_cox$glmnet_model <- glmnet(x = multi_cox$x$tcga, 
                                   y = multi_cox$y$tcga, 
                                   family = 'cox', 
                                   alpha = multi_cox$alpha, 
                                   lambda = multi_cox$opt_lambda$lambda)
  
# Linear predictor scores ------
  
  insert_msg('Linear predictor scores')
  
  multi_cox$lp_scores <- multi_cox$x %>% 
    map(~predict(multi_cox$glmnet_model, 
                 newx = .x, 
                 type = 'link')) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'patient_id') %>% 
    map(set_names, c('patient_id', 'lp_score')) %>% 
    map(as_tibble) %>% 
    map2(., multi_cox$analysis_tbl, 
         left_join, by = 'patient_id')
  
  ## stratification of the lp score in tertiles 
  
  multi_cox$lp_scores <- multi_cox$lp_scores %>%
    map(~mutate(.x, 
                lp_strata = cut(lp_score, 
                                c(-Inf, quantile(lp_score, c(.33, 0.66)), Inf), 
                                c('low', 'intermediate', 'high'))))
  
# Uni variate Cox models of the LP score ------
  
  insert_msg('Uni-variate Cox models')
  
  multi_cox$cox_models <- multi_cox$lp_scores %>% 
    map(~coxph(Surv(os_days, death) ~ lp_score, 
               data = .x)) %>% 
    map2(., multi_cox$analysis_tbl, 
         as_coxex)
  
# Assumptions ------
  
  insert_msg('Assumptions')
  
  multi_cox$assumptions <- multi_cox$cox_models %>% 
    map(summary, 'assumptions')
  
# Fit stats -------
  
  insert_msg('Fit stats')
  
  multi_cox$fit_stats <- multi_cox$cox_models %>% 
    map(summary, 'fit') %>%
    compress(names_to = 'cohort') %>% 
    mutate(plot_cap = paste0('R\u00B2 = ', signif(raw_rsq, 2), 
                             ', C = ', signif(c_index, 2), 
                             ' [', signif(lower_ci, 2), ' - ', 
                             signif(upper_ci, 2), ']'), 
           plot_tag = paste0('\ntotal: n = ', n_complete, 
                            ', events: n = ', n_events))
  
# Non-zero model coefficients ------
  
  insert_msg('Non-zero coefficients')
  
  multi_cox$coefs <- multi_cox$glmnet_model %>% 
    coef %>% 
    as.matrix %>% 
    as.data.frame %>% 
    rownames_to_column('variable') %>% 
    set_names(c('variable', 'hr')) %>% 
    filter(hr != 0) %>% 
    mutate(hr = exp(hr))
  
# Coefficient plot ---------
  
  insert_msg('Bubble plot with the coefficients')
  
  multi_cox$coef_plot <- multi_cox$coefs %>% 
    ggplot(aes(x = log(hr), 
               y = reorder(variable, hr), 
               size = abs(log(hr)), 
               fill = factor(sign(log(hr))))) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_point(shape = 21) + 
    geom_text(aes(label = signif(hr, 3)), 
              size = 2.6, 
              vjust = -1.6, 
              hjust= 0.5) + 
    scale_fill_manual(values = c('-1' = 'steelblue', 
                                 '0' = 'gray60', 
                                 '1' = 'firebrick')) + 
    guides(fill = 'none') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(face = 'italic')) + 
    labs(title = 'Non-zero model coefficients, TCGA KIRC', 
         subtitle = 'Training cohort', 
         x = expression('HR'[Ridge]))
  
# LP score tertile KM plots --------- 
  
  insert_msg('Tertile plots of the coefficients')
  
  ## survfit objects
  
  multi_cox$surv_fit <- multi_cox$lp_scores %>% 
    map(~surv_fit(Surv(os_days, death) ~ lp_strata, 
                  data = .x))
  
  ## log-rank test
  
  multi_cox$surv_diff <- multi_cox$surv_fit %>% 
    map(surv_pvalue, 
        method = 'survdiff', 
        test.for.trend = TRUE) %>% 
    map(as_tibble) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(p_adjusted = p.adjust(pval, 'BH'), 
           plot_lab = ifelse(p_adjusted < 0.05, 
                             paste('p =', signif(p_adjusted, 2)),
                             paste0('ns (p = ', signif(p_adjusted, 2), ')')))
  
  ## strata n numbers
  
  multi_cox$surv_strata_n <- multi_cox$lp_scores %>% 
    map(count, lp_strata) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ', n = '))
  
  ## plots
  
  multi_cox$km_plots <- 
    list(fit = multi_cox$surv_fit, 
         title = globals$cohort_lexicon$label[names(multi_cox$surv_fit)], 
         legend.labs = multi_cox$surv_strata_n, 
         pval = multi_cox$surv_diff$plot_lab) %>% 
    pmap(ggsurvplot, 
         palette = c('darkolivegreen4', 
                     'steelblue3', 
                     'coral3'), 
         legend.title = 'Ridge OS score', 
         xlab = 'Overall survival, days', 
         pval.size = 2.75) %>% 
    map(~.x$plot + 
          globals$common_theme)
  
  ## appending with the fit statistics (R-square, C-index)
  ## and numbers of cases
  
  multi_cox$km_plots <- 
    list(x = multi_cox$km_plots, 
         y = multi_cox$fit_stats$plot_cap, 
         z = multi_cox$fit_stats$plot_tag) %>% 
    pmap(function(x, y, z) x + 
           labs(subtitle = y, 
                tag = z))
  
# END ------
  
  plan('sequential')
  
  rm(i)
  
  insert_tail()