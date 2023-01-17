# Comparison of overall survival and relapse-free survival between the clusters.

  insert_head()
  
# container ---------
  
  clust_surv <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## analysis tables
  
  clust_surv$analysis_tbl <- list()
  
  for(i in globals$os_cohorts) {
    
    eval(parse_expr(paste0('clust_surv$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
  }
  
  clust_surv$analysis_tbl$tcga <- clust_surv$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  clust_surv$analysis_tbl <- clust_surv$analysis_tbl %>% 
    map(select, 
        patient_id, 
        any_of(globals$surv_vars))
  
  ## appending with the cluster assignment
  
  clust_surv$analysis_tbl <- gene_clust$clust_obj[globals$os_cohorts] %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(clust_surv$analysis_tbl, ., 
         left_join, by = 'patient_id')
  
# N numbers ------
  
  insert_msg('N numbers')
  
  ## OS
  
  clust_surv$n_numbers$os <- clust_surv$analysis_tbl %>% 
    map(select, clust_id, os_days, death) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(count, clust_id)
  
  clust_surv$plot_labs$os <- clust_surv$n_numbers$os %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = '))
  
  ## RFS
  
  clust_surv$n_numbers$rfs <- 
    clust_surv$analysis_tbl[c('tcga', 'cm10', 
                              'cm25ev', 'cm25ni')] %>% 
    map(select, clust_id, rfs_days, relapse) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(count, clust_id)
  
  clust_surv$plot_labs$rfs <- clust_surv$n_numbers$rfs %>% 
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = '))
  
# Numbers of total observations and the events ------
  
  insert_msg('Numbers of observation and events')
  
  ## OS
  
  clust_surv$n_events$os <- clust_surv$analysis_tbl %>% 
    map(select, clust_id, os_days, death) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(count, death)
  
  clust_surv$plot_caps$os <- clust_surv$n_events$os %>% 
    map(~paste0('total: n = ', 
                sum(.x[[2]]), 
                ', events: n = ',
                .x[[2]][2]))
  
  ## RFS
  
  clust_surv$n_events$rfs <- 
    clust_surv$analysis_tbl[c('tcga', 'cm10', 
                              'cm25ev', 'cm25ni')] %>% 
    map(select, clust_id, rfs_days, relapse) %>% 
    map(~filter(.x, complete.cases(.x))) %>% 
    map(count, relapse)
  
  clust_surv$plot_caps$rfs <- clust_surv$n_events$rfs %>% 
    map(~paste0('total: n = ', 
                sum(.x[[2]]), 
                ', events: n = ',
                .x[[2]][2]))
  
# Survfit objects --------
  
  insert_msg('Sur fit objects')
  
  clust_surv$survfit_obj$os <- clust_surv$analysis_tbl %>% 
    map(~surv_fit(Surv(os_days, death) ~ clust_id, 
                  data = .x))
  
  clust_surv$survfit_obj$rfs <- 
    clust_surv$analysis_tbl[c('tcga', 
                              'cm10', 
                              'cm25ev', 
                              'cm25ni')] %>% 
    map(~surv_fit(Surv(rfs_days, relapse) ~ clust_id, 
                  data = .x))
  
# Significance of the survival differences ------
  
  insert_msg('Significance of the survival differences')
  
  clust_surv$test_results <- clust_surv$survfit_obj  %>% 
    map(~map_dfr(.x, 
                 surv_pvalue, 
                 method = 'S2')) %>% 
    map(mutate, 
        p_adjusted = p.adjust(pval, 'BH'), 
        sign_raw = ifelse(p_adjusted < 0.05, 
                          paste('raw: p =', signif(pval, 2)), 
                          paste0('raw: ns (p = ', signif(pval, 2), ')')), 
        sign_fdr = ifelse(p_adjusted < 0.05, 
                          paste('FDR: p =', signif(p_adjusted, 2)), 
                          paste0('FDR: ns (p = ', signif(p_adjusted, 2), ')')), 
        significance = paste(sign_raw, sign_fdr, sep = '\n'))
  
# Kaplan-Meier plots -------
  
  insert_msg('Kaplan Meier plots')
  
  ## overall survival
  
  clust_surv$km_plots$os <- 
    list(fit = clust_surv$survfit_obj$os,
         pval = clust_surv$test_results$os$significance, 
         title = globals$cohort_lexicon$label[names(clust_surv$survfit_obj$os)], 
         legend.labs = clust_surv$plot_labs$os) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors), 
         xlab = 'Overall survival, days', 
         legend.title = 'Lymphokine cluster', 
         pval.size = 2.75) %>% 
    map2(., clust_surv$plot_caps$os, 
         ~.x$plot + 
           labs(subtitle = .y) + 
           globals$common_theme)
  
  ## relapse-free survival
  
  clust_surv$km_plots$rfs <- 
    list(fit = clust_surv$survfit_obj$rfs,
         pval = clust_surv$test_results$rfs$significance, 
         title = globals$cohort_lexicon$label[names(clust_surv$survfit_obj$rfs)], 
         legend.labs = clust_surv$plot_labs$rfs) %>% 
    pmap(ggsurvplot, 
         palette = unname(globals$cluster_colors), 
         xlab = 'Relapse-free survival, days', 
         legend.title = 'Lymphokine cluster', 
         pval.size = 2.75) %>% 
    map2(., clust_surv$plot_caps$rfs, 
         ~.x$plot + 
           labs(subtitle = .y) + 
           globals$common_theme)

# END ------
  
  rm(i)