# Multi-parameter modeling of overall survival in the three major cohorts:
# TCGA KIRC, E-MTAB 1980 amd RECA-EU with the genes of interest 
# and base clinical paramaters (age, sex, grade, staging)
#
# The TCGA cohort serves as a training collective. The linear predictor scores
# of the TCGA GLMNET model are predicted for the remaining collectives 
# and the association of the linear predictor scores with survival is modeled
# by uni-variable Cox regression


  insert_head()
  
# container -----
  
  clin_os <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## variables: demographic and clinical features 
  ## present in all cohorts, genes of interest
  ## node stage is not included since it is highly missing in the TCGA cohort
  
  clin_os$genes <- globals$cxcl_genes
  
  clin_os$variables <- c('age', 'sex', 
                         'tumor_grade', 
                         'pt_stage', 'pm_stage', 
                         clin_os$genes)
  
  ## analysis tables with the complete set of variables
  ## recoding of the staging and grades (e.g. TX turned into NA, simplification
  ## of the tumor grades)
  
  clin_os$analysis_tbl <- list()
  
  for(i in c('tcga', 'emtab1980', 'reca')) {
    
    eval(parse_expr(paste0('clin_os$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  clin_os$analysis_tbl <- clin_os$analysis_tbl %>% 
    map(select, 
        patient_id, 
        os_days,
        death, 
        any_of('tissue_type'), 
        all_of(clin_os$variables))
  
  clin_os$analysis_tbl$tcga <- clin_os$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  clin_os$analysis_tbl <- clin_os$analysis_tbl %>% 
    map(mutate, 
        tumor_grade = stri_extract(tumor_grade, regex = '\\d{1}'), 
        tumor_grade = paste0('G', tumor_grade), 
        tumor_grade = ifelse(tumor_grade %in% c('G1', 'G2'), 
                             'G1 - G2', 
                             ifelse(tumor_grade %in% c('G3', 'G4'), 
                                    'G3 - G4', NA)), 
        tumor_grade = factor(tumor_grade, c('G1 - G2', 'G3 - G4')), 
        pt_stage = stri_extract(pt_stage, regex = 'T\\d{1}'), 
        pt_stage = ifelse(pt_stage %in% c('T1', 'T2'), 
                          'T1 - T2', 
                          ifelse(pt_stage %in% c('T3', 'T4'), 
                                 'T3 - T4', NA)), 
        pt_stage = factor(pt_stage, c('T1 - T2', 'T3 - T4')), 
        pm_stage = factor(pm_stage, c('M0', 'M1')))
  
  clin_os$analysis_tbl <- clin_os$analysis_tbl %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## adding a day to the OS, immediate censoring is not allowed for GLMNET
  
  clin_os$analysis_tbl <- clin_os$analysis_tbl %>% 
    map(mutate, os_days = os_days + 1)
  
  ## normalization of the numeric explanatory features
  ## second order term for age
  
  for(i in c('age', clin_os$genes)) {
    
    clin_os$analysis_tbl <- clin_os$analysis_tbl %>% 
      map(mutate, 
          !!i := scale(.data[[i]])[, 1], 
          age_sq = age^2)
    
  }
  
  clin_os$variables <- c(clin_os$variables, 'age_sq')
  
  ## survival objects and explanatory variable matrices
  ## no higher order terms for gene expression variables 
  ## included in the models
  ## since no significant effects for splines were detected
  ## in uni-variable Cox modeling
  ##
  ## creating a modeling matrix without gene variables. 
  ## It will be used to check the added effect of gene expression
  ## to the accuracy of survival prediction
  
  clin_os$y <- clin_os$analysis_tbl %>% 
    map(~Surv(.x$os_days, .x$death))
  
  clin_os$x <- clin_os$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(~.x[clin_os$variables]) %>% 
    map(~model.matrix(~., data = .x))
  
  clin_os$x_clin <- clin_os$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(~.x[clin_os$variables[!clin_os$variable %in% clin_os$genes]]) %>% 
    map(~model.matrix(~., data = .x))
  
  ## alpha and CV folds
  
  clin_os$alpha <- 0 
  clin_os$n_rep <- 200
  
  ## folds for CV: TCGA is the training cohort
  
  set.seed(1234)
  
  clin_os$folds <- 1:clin_os$n_rep %>% 
    map(~createFolds(y = clin_os$analysis_tbl$tcga$patient_id, 
                     k = 10, 
                     list = FALSE, 
                     returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:clin_os$n_rep))

# lambda tuning -------
  
  insert_msg('Lambda tuning, TCGA')
  
  ## model construction: with the genes
  
  clin_os$lambda_tune$genes <- clin_os$folds %>% 
    future_map(~cv.glmnet(x = clin_os$x$tcga, 
                          y = clin_os$y$tcga, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = clin_os$alpha), 
               .options = furrr_options(seed = TRUE))
  
  clin_os$lambda_tune$clin <- clin_os$folds %>% 
    future_map(~cv.glmnet(x = clin_os$x_clin$tcga, 
                          y = clin_os$y$tcga, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = clin_os$alpha), 
               .options = furrr_options(seed = TRUE))
  
  ## extracting the optimal lambdas
  
  clin_os$lambda_tbl <- clin_os$lambda_tune %>% 
    map(~map(.x, 
             ~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo', 
                             'lambda.min', 'lambda.1se')])) %>% 
          map_dfr(filter, lambda == lambda.1se))
  
  clin_os$opt_lambda <- clin_os$lambda_tbl %>% 
    map(filter, cvm == min(cvm))
  
# Building the GLMNET models -------
  
  insert_msg('Building the GLMNET models')
  
  ## the full model including gene expression variables
  ## and the clinical variable-only model
  
  clin_os$glmnet_model$genes <- 
    glmnet(x = clin_os$x$tcga, 
           y = clin_os$y$tcga, 
           family = 'cox', 
           alpha = clin_os$alpha, 
           lambda = clin_os$opt_lambda$genes$lambda)
  
  clin_os$glmnet_model$clin <- 
    glmnet(x = clin_os$x_clin$tcga, 
           y = clin_os$y$tcga, 
           family = 'cox', 
           alpha = clin_os$alpha, 
           lambda = clin_os$opt_lambda$clin$lambda)
  
# Linear predictor scores ------
  
  insert_msg('Linear predictor scores')
  
  ## predictions: the full model and the clinical variable-only model
  
  clin_os$lp_scores$genes <- clin_os$x %>% 
    map(~predict(clin_os$glmnet_model$genes, 
                 newx = .x, 
                 type = 'link'))
  
  clin_os$lp_scores$clin <- clin_os$x_clin %>% 
    map(~predict(clin_os$glmnet_model$clin, 
                 newx = .x, 
                 type = 'link'))
  
  ## appending with patient ID and expression information
  
  clin_os$lp_scores <- clin_os$lp_scores %>% 
    map(~map(.x, as.data.frame) %>% 
          map(rownames_to_column, 'patient_id') %>% 
          map(set_names, c('patient_id', 'lp_score')) %>% 
          map(as_tibble) %>% 
          map2(., clin_os$analysis_tbl, 
               left_join, by = 'patient_id'))

  ## stratification of the lp score in tertiles 
  
  clin_os$lp_scores <- clin_os$lp_scores %>%
    map(~map(.x, 
             ~mutate(.x, 
                lp_strata = cut(lp_score, 
                                c(-Inf, quantile(lp_score, c(.33, 0.66)), Inf), 
                                c('low', 'intermediate', 'high')))))
  
# Uni-variate Cox models of the LP score ------
  
  insert_msg('Uni-variate Cox models')
  
  ## for the full models and the clinical variable-only models
  
  clin_os$cox_models$genes <- clin_os$lp_scores$genes %>% 
    map(~coxph(Surv(os_days, death) ~ lp_score, 
               data = .x)) %>% 
    map2(., clin_os$analysis_tbl, 
         as_coxex)
  
  clin_os$cox_models$clin <- clin_os$lp_scores$clin %>% 
    map(~coxph(Surv(os_days, death) ~ lp_score, 
               data = .x)) %>% 
    map2(., clin_os$analysis_tbl, 
         as_coxex)
  
# Assumptions ------
  
  insert_msg('Assumptions')
  
  clin_os$assumptions <- clin_os$cox_models %>% 
    map(~map(.x, summary, 'assumptions'))
  
# Fit stats -------
  
  insert_msg('Fit stats')
  
  clin_os$fit_stats <- clin_os$cox_models %>% 
    map(~map(.x, summary, 'fit') %>%
          compress(names_to = 'cohort') %>% 
          mutate(plot_cap = paste0('R\u00B2 = ', signif(raw_rsq, 2), 
                                   ', C = ', signif(c_index, 2), 
                                   ' [', signif(lower_ci, 2), ' - ', 
                                   signif(upper_ci, 2), ']'), 
                 plot_tag = paste0('\ntotal: n = ', n_complete, 
                                   ', events: n = ', n_events), 
                 cohort_lab = exchange(cohort, 
                                       dict = globals$cohort_lexicon, 
                                       key = 'cohort', 
                                       value = 'label'), 
                 ax_lab = paste0(cohort_lab, 
                                 '\ntotal: n = ', n_complete, 
                                 '\nevents: n = ', n_events)))
  
# Fit stat plots for the full models and clinical-only models -----
  
  insert_msg('Fit stat plots')
  
  ## C-index
  
  clin_os$c_index_plot <- clin_os$fit_stats %>% 
    compress(names_to = 'model_type') %>% 
    ggplot(aes(x = c_index, 
               y = cohort, 
               color = model_type)) +
    geom_vline(xintercept = 0.5, linetype = 'dashed') + 
    geom_errorbarh(aes(xmin = lower_ci, 
                       xmax = upper_ci),
                   height = 0, 
                   position = position_dodge(0.75)) + 
    geom_point(shape = 16, 
               position = position_dodge(0.75)) + 
    geom_text(aes(label = paste0(signif(c_index, 2), 
                                 ' [', signif(lower_ci, 2), 
                                 ' - ', signif(upper_ci, 2), ']')), 
              size = 2.6, 
              hjust = 0.5, 
              vjust = -1.3,
              position = position_dodge(0.75)) + 
    scale_x_continuous(limits = c(0.5, 1)) + 
    scale_color_manual(values = c(clin = 'steelblue3', 
                                  genes = 'coral3'), 
                       labels = c(clin = 'clinical-only', 
                                  genes = 'genes + clinical'), 
                       name = 'Model type') + 
    scale_y_discrete(labels =  set_names(clin_os$fit_stats$genes$ax_lab, 
                                         clin_os$fit_stats$genes$cohort)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Goodness of fit, multi-paramater modeling', 
         subtitle = 'Ridge Cox regression', 
         x = 'C-index, 95% CI')
  
# Non-zero model coefficients ------
  
  insert_msg('Non-zero coefficients')
  
  clin_os$coefs <- clin_os$glmnet_model %>% 
    map(coef) %>% 
    map(as.matrix) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'variable') %>% 
    map(set_names, c('parameter', 'hr')) %>% 
    map(filter, hr != 0) %>% 
    map(mutate, hr = exp(hr))
  
  ## appending with the variable labels and levels
  
  clin_os$coefs <- clin_os$coefs %>% 
    map(mutate, 
        variable = ifelse(parameter == 'age_sq', 
                          'age_sq', 
                          stri_extract(parameter, 
                                       regex = paste(clin_os$variables, 
                                                     collapse = '|'))), 
        level = ifelse(parameter == 'age_sq', 
                       NA, 
                       stri_replace(parameter, 
                                    regex = paste(clin_os$variables, 
                                                  collapse = '|'), 
                                    replacement = '')), 
        level = ifelse(level == '', NA, level), 
        var_lab = ifelse(variable == 'age_sq', 
                         'Age\u00B2', 
                         exchange(variable, 
                                  dict = globals$var_lexicon, 
                                  key = 'variable', 
                                  value = 'label')), 
        var_lab = ifelse(variable %in% clin_os$genes, 
                         paste0('<em>', variable, '</em>'), 
                         var_lab), 
        ax_lab = ifelse(is.na(level), 
                        var_lab, 
                        paste(var_lab, level, sep = ': ')))
  
# Coefficient plot ---------
  
  insert_msg('Bubble plot with the coefficients')
  
  clin_os$coef_plots <- 
    map2(clin_os$coefs, 
         c('Clinical variable and gene model', 
           'Clnical variable-only model'), 
         ~ggplot(.x, 
                 aes(x = log(hr), 
                     y = reorder(ax_lab, hr), 
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
                 axis.text.y = element_markdown()) + 
           labs(title = 'Model coefficients, TCGA KIRC', 
                subtitle = .y, 
                x = expression('log HR'[Ridge])))
    
# LP score tertile KM plots: for the full models only --------- 
  
  insert_msg('Tertile plots of the coefficients')
  
  ## survfit objects
  
  clin_os$surv_fit <- clin_os$lp_scores$genes %>% 
    map(~surv_fit(Surv(os_days, death) ~ lp_strata, 
                  data = .x))
  
  ## log-rank test
  
  clin_os$surv_diff <- clin_os$surv_fit %>% 
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
  
  clin_os$surv_strata_n <- clin_os$lp_scores$genes %>% 
    map(count, lp_strata) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ', n = '))
  
  ## plots
  
  clin_os$km_plots <- 
    list(fit = clin_os$surv_fit, 
         title = globals$cohort_lexicon$label[names(clin_os$surv_fit)], 
         legend.labs = clin_os$surv_strata_n, 
         pval = clin_os$surv_diff$plot_lab) %>% 
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
  
  clin_os$km_plots <- 
    list(x = clin_os$km_plots, 
         y = clin_os$fit_stats$genes$plot_cap, 
         z = clin_os$fit_stats$genes$plot_tag) %>% 
    pmap(function(x, y, z) x + 
           labs(subtitle = y, 
                tag = z))
    
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()