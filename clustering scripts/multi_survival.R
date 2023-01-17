# Multi-parameter modeling of overall survival in the lymphokine clusters
# Clinical and demographic parameters (age, sex, grade and stage) serve along
# with the cluster assignment information as explanatory variables
#
# Modeling is done with Ridge Cox regression

  insert_head()
  
# container ------
  
  multi_clust <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## variables. Node staging is excluded, bicause its highly missing 
  ## in the TCGA training cohort
  
  multi_clust$variables <- c('age', 'sex', 
                             'tumor_grade', 
                             'pt_stage', 'pm_stage')
  
  ## analysis tables with the complete set of variables
  ## recoding of the staging and grades (e.g. TX turned into NA, simplification
  ## of the tumor grades)
  ## normaliztion of age and second order
  
  multi_clust$analysis_tbl <- list()
  
  for(i in c('tcga', 'emtab1980', 'reca')) {
    
    eval(parse_expr(paste0('multi_clust$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }

  multi_clust$analysis_tbl$tcga <- multi_clust$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  
  multi_clust$analysis_tbl <- multi_clust$analysis_tbl %>% 
    map(select, 
        patient_id, 
        os_days,
        death, 
        all_of(multi_clust$variables))
  
  multi_clust$analysis_tbl <- multi_clust$analysis_tbl %>% 
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
        pm_stage = factor(pm_stage, c('M0', 'M1')), 
        age = scale(age)[, 1], 
        age_sq = age^2)
  
  multi_clust$variables <- c(multi_clust$variables, 'age_sq')

  ## adding a day to the OS, immediate censoring is not allowed for GLMNET
  
  multi_clust$analysis_tbl <- multi_clust$analysis_tbl %>% 
    map(mutate, os_days = os_days + 1)
  
  ## appending the analysis tables with the cluster assignment scheme
  ## complete cases
  
  multi_clust$analysis_tbl <- 
    gene_clust$clust_obj[names(multi_clust$analysis_tbl)] %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(., multi_clust$analysis_tbl, 
         left_join, by = 'patient_id')
  
  multi_clust$analysis_tbl <- multi_clust$analysis_tbl %>% 
    map(~filter(.x, complete.cases(.x)))
  
  ## survival objects and matrices of explanatory variables
  ##
  ## creating a modeling matrix without gene variables. 
  ## It will be used to check the added effect of the clustering
  ## to the accuracy of survival prediction
  
  multi_clust$y <- multi_clust$analysis_tbl %>% 
    map(~Surv(.x$os_days, .x$death))
  
  multi_clust$x <- multi_clust$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(~.x[c(multi_clust$variables, 'clust_id')]) %>% 
    map(~model.matrix(~., data = .x))
  
  multi_clust$x_clin <- multi_clust$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(~.x[multi_clust$variables]) %>% 
    map(~model.matrix(~., data = .x))
  
  multi_clust$variables <- c(multi_clust$variables, 'clust_id')
  
  ## alpha and CV folds
  
  multi_clust$alpha <- 0 
  multi_clust$n_rep <- 200
  
  ## folds for CV: TCGA is the training cohort
  
  set.seed(1234)
  
  multi_clust$folds <- 1:multi_clust$n_rep %>% 
    map(~createFolds(y = multi_clust$analysis_tbl$tcga$patient_id, 
                     k = 10, 
                     list = FALSE, 
                     returnTrain = TRUE)) %>% 
    set_names(paste0('rep_', 1:multi_clust$n_rep))
  
# lambda tuning -------
  
  insert_msg('Lambda tuning, TCGA')
  
  ## model construction: with the genes
  
  multi_clust$lambda_tune$cluster <- multi_clust$folds %>% 
    future_map(~cv.glmnet(x = multi_clust$x$tcga, 
                          y = multi_clust$y$tcga, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = multi_clust$alpha), 
               .options = furrr_options(seed = TRUE))
  
  multi_clust$lambda_tune$clin <- multi_clust$folds %>% 
    future_map(~cv.glmnet(x = multi_clust$x_clin$tcga, 
                          y = multi_clust$y$tcga, 
                          foldid = .x, 
                          type.measure = 'default', 
                          family = 'cox', 
                          alpha = multi_clust$alpha), 
               .options = furrr_options(seed = TRUE))
  
  ## extracting the optimal lambdas
  
  multi_clust$lambda_tbl <- multi_clust$lambda_tune %>% 
    map(~map(.x, 
             ~as_tibble(.x[c('lambda', 'cvm', 'cvup', 'cvlo', 
                             'lambda.min', 'lambda.1se')])) %>% 
          map_dfr(filter, lambda == lambda.1se))
  
  multi_clust$opt_lambda <- multi_clust$lambda_tbl %>% 
    map(filter, cvm == min(cvm))
  
# Building the GLMNET models -------
  
  insert_msg('Building the GLMNET models')
  
  ## the full model including clustering
  ## and the clinical variable-only model
  
  multi_clust$glmnet_model$cluster <- 
    glmnet(x = multi_clust$x$tcga, 
           y = multi_clust$y$tcga, 
           family = 'cox', 
           alpha = multi_clust$alpha, 
           lambda = multi_clust$opt_lambda$cluster$lambda)
  
  multi_clust$glmnet_model$clin <- 
    glmnet(x = multi_clust$x_clin$tcga, 
           y = multi_clust$y$tcga, 
           family = 'cox', 
           alpha = multi_clust$alpha, 
           lambda = multi_clust$opt_lambda$clin$lambda)
  
# Linear predictor scores ------
  
  insert_msg('Linear predictor scores')
  
  ## predictions: the full model and the clinical variable-only model
  
  multi_clust$lp_scores$cluster <- multi_clust$x %>% 
    map(~predict(multi_clust$glmnet_model$cluster, 
                 newx = .x, 
                 type = 'link'))
  
  multi_clust$lp_scores$clin <- multi_clust$x_clin %>% 
    map(~predict(multi_clust$glmnet_model$clin, 
                 newx = .x, 
                 type = 'link'))
  
  ## appending with patient ID and expression information
  
  multi_clust$lp_scores <- multi_clust$lp_scores %>% 
    map(~map(.x, as.data.frame) %>% 
          map(rownames_to_column, 'patient_id') %>% 
          map(set_names, c('patient_id', 'lp_score')) %>% 
          map(as_tibble) %>% 
          map2(., multi_clust$analysis_tbl, 
               left_join, by = 'patient_id'))
  
  ## stratification of the lp score in tertiles 
  
  multi_clust$lp_scores <- multi_clust$lp_scores %>%
    map(~map(.x, 
             ~mutate(.x, 
                     lp_strata = cut(lp_score, 
                                     c(-Inf, 
                                       quantile(lp_score, c(.33, 0.66)), 
                                       Inf), 
                                     c('low', 
                                       'intermediate', 
                                       'high')))))
  
# Uni-variate Cox models of the LP score ------
  
  insert_msg('Uni-variate Cox models')
  
  ## for the full models and the clinical variable-only models
  
  multi_clust$cox_models$cluster <- multi_clust$lp_scores$cluster %>% 
    map(~coxph(Surv(os_days, death) ~ lp_score, 
               data = .x)) %>% 
    map2(., multi_clust$analysis_tbl, 
         as_coxex)
  
  multi_clust$cox_models$clin <- multi_clust$lp_scores$clin %>% 
    map(~coxph(Surv(os_days, death) ~ lp_score, 
               data = .x)) %>% 
    map2(., multi_clust$analysis_tbl, 
         as_coxex)
  
# Assumptions ------
  
  insert_msg('Assumptions')
  
  multi_clust$assumptions <- multi_clust$cox_models %>% 
    map(~map(.x, summary, 'assumptions'))
  
# Fit stats -------
  
  insert_msg('Fit stats')
  
  multi_clust$fit_stats <- multi_clust$cox_models %>% 
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
                                 ', events: n = ', n_events)))
  
# Fit stat plots for the full models and clinical-only models -----
  
  insert_msg('Fit stat plots')
  
  ## C-index
  
  multi_clust$c_index_plot <- multi_clust$fit_stats %>% 
    compress(names_to = 'model_type') %>% 
    ggplot(aes(x = c_index, 
               y = cohort, 
               color = model_type)) + 
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
              vjust = -1.5,
              position = position_dodge(0.75)) + 
    scale_color_manual(values = c(clin = 'steelblue3', 
                                  cluster = 'coral3'), 
                       labels = c(clin = 'clinical-only', 
                                  luster = 'lymphokine clusters + clinical'), 
                       name = 'Model type') + 
    scale_y_discrete(labels =  set_names(multi_clust$fit_stats$cluster$ax_lab, 
                                         multi_clust$fit_stats$cluster$cohort)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Goodness of fit, multi-paramater modeling', 
         subtitle = 'Ridge Cox regression', 
         x = 'C-index, 95% CI')
  
# Non-zero model coefficients ------
  
  insert_msg('Non-zero coefficients')
  
  multi_clust$coefs <- multi_clust$glmnet_model %>% 
    map(coef) %>% 
    map(as.matrix) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'variable') %>% 
    map(set_names, c('parameter', 'hr')) %>% 
    map(filter, hr != 0) %>% 
    map(mutate, hr = exp(hr))
  
  ## appending with the variable labels and levels
  
  multi_clust$coefs <- multi_clust$coefs %>% 
    map(mutate, 
        variable = ifelse(parameter == 'age_sq', 
                          'age_sq', 
                          stri_extract(parameter, 
                                       regex = paste(multi_clust$variables, 
                                                     collapse = '|'))), 
        level = ifelse(parameter == 'age_sq', 
                       NA, 
                       stri_replace(parameter, 
                                    regex = paste(multi_clust$variables, 
                                                  collapse = '|'), 
                                    replacement = '')), 
        level = ifelse(level == '', NA, level), 
        var_lab = ifelse(variable == 'age_sq', 
                         'Age\u00B2', 
                         exchange(variable, 
                                  dict = globals$var_lexicon, 
                                  key = 'variable', 
                                  value = 'label')), 
        var_lab = ifelse(variable == 'clust_id', 
                         'Lymphokine cluster', 
                         var_lab), 
        ax_lab = ifelse(is.na(level), 
                        var_lab, 
                        paste(var_lab, level, sep = ': ')))
  
# Coefficient plot ---------
  
  insert_msg('Bubble plot with the coefficients')
  
  multi_clust$coef_plots <- 
    map2(multi_clust$coefs, 
         c('Clinical variable and lymphokine cluster model', 
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
           labs(title = 'Non-zero model coefficients, TCGA KIRC', 
                subtitle = .y, 
                x = expression('HR'[Ridge])))
  
# LP score tertile KM plots: for the full models only --------- 
  
  insert_msg('Tertile plots of the coefficients')
  
  ## survfit objects
  
  multi_clust$surv_fit <- multi_clust$lp_scores$cluster %>% 
    map(~surv_fit(Surv(os_days, death) ~ lp_strata, 
                  data = .x))
  
  ## log-rank test
  
  multi_clust$surv_diff <- multi_clust$surv_fit %>% 
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
  
  multi_clust$surv_strata_n <- multi_clust$lp_scores$cluster %>% 
    map(count, lp_strata) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ', n = '))
  
  ## plots
  
  multi_clust$km_plots <- 
    list(fit = multi_clust$surv_fit, 
         title = globals$cohort_lexicon$label[names(multi_clust$surv_fit)], 
         legend.labs = multi_clust$surv_strata_n, 
         pval = multi_clust$surv_diff$plot_lab) %>% 
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
  
  multi_clust$km_plots <- 
    list(x = multi_clust$km_plots, 
         y = multi_clust$fit_stats$cluster$plot_cap, 
         z = multi_clust$fit_stats$cluster$plot_tag) %>% 
    pmap(function(x, y, z) x + 
           labs(subtitle = y, 
                tag = z))
  
# END -----
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()