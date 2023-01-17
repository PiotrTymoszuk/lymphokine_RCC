# Analysis of overall survival in low-stage (T1 - T2) and  high-stage cancers
# stratified in low and high expressors of chemokine genes
#
# Done for the TCGA, e-MTAB-1980 and RECA-EU collectives

  insert_head()
  
# container --------
  
  stage_cut <- list()
  
# parallel backend --------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## analysis cohorts
  
  stage_cut$cohorts <- c('tcga', 'emtab1980', 'reca')
  
  stage_cut$cohort_labs <- stage_cut$cohorts %>% 
    exchange(dict = globals$cohort_lexicon, 
             key = 'cohort',
             value = 'label') %>% 
    map(~paste(.x, c('T1 - T2', 'T3 - T4'))) %>% 
    unlist
  
  ## analysis tables with normalized expression
  ## of the genes of interest
  
  stage_cut$analysis_tbl <- list()
  
  for(i in stage_cut$cohorts) {
    
    eval(parse_expr(paste0('stage_cut$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  stage_cut$analysis_tbl <- stage_cut$analysis_tbl %>% 
    map(select, 
        patient_id, 
        pt_stage, 
        os_days,
        death, 
        any_of('tissue_type'), 
        all_of(globals$cxcl_genes))
  
  stage_cut$analysis_tbl$tcga <- stage_cut$analysis_tbl$tcga %>% 
    filter(tissue_type == 'Tumor')
  
  for(i in globals$cxcl_genes) {
    
    stage_cut$analysis_tbl <- stage_cut$analysis_tbl %>% 
      map(mutate, 
          !!i := scale(.data[[i]])[, 1])
    
  }
  
  ## splitting by the stage
  
  stage_cut$analysis_tbl <- stage_cut$analysis_tbl %>% 
    map(mutate, 
        pt_stage = stri_extract(pt_stage, regex = 'T\\d{1}'), 
        progression = ifelse(pt_stage %in% c('T1', 'T2'), 
                             'early', 'late'), 
        progression = factor(progression, c('early', 'late'))) %>% 
    map(~filter(.x, complete.cases(.x)))
  
  stage_cut$analysis_tbl <- stage_cut$analysis_tbl %>% 
    map(dlply, 'progression', as_tibble) %>% 
    unlist(recursive = FALSE)
  
  ## minimal strata size (20% of the all observations)
  
  stage_cut$min_n <- stage_cut$analysis_tbl %>% 
    map_dbl(nrow)
  
  stage_cut$min_n <- 0.2 * stage_cut$min_n

# finding the cutoffs -------
  
  insert_msg('Stratification')
  
  stage_cut$cut_obj[globals$cxcl_genes] <- 
    globals$cxcl_genes %>%
    map(function(gene) list(data = stage_cut$analysis_tbl,
                            min_n = stage_cut$min_n) %>%
          future_pmap(find_cutoff,
                      time = 'os_days', 
                      event = 'death', 
                      variable = gene, 
                      rho = 0,
                      .parallel = FALSE))
  
# Optimal cutoff stats, p value correction within the cohort ------
  
  insert_msg('Cutoff summary, p correction')
  
  stage_cut$test_results <-
    stage_cut$cut_obj %>%
    map(~map2_dfr(.x, names(.x),
                  ~mutate(summary(.x)[1, ], model = .y)))
  
  ## stripping the information to be shown in the plot captions
  ## gene-wise correction form multiple testing
  
  stage_cut$test_results <- stage_cut$test_results %>%
    map(mutate,
        p_adjusted = p.adjust(p_value, 'BH'),
        plot_cap = paste0('pFDR = ', signif(p_adjusted, 2),
                          ', p raw = ', signif(p_value, 2)),
        plot_cap = ifelse(p_adjusted < 0.05,
                          plot_cap,
                          paste('ns,', plot_cap)),
        plot_cap = paste0('cutoff = ',
                          signif(cutoff, 2),
                          ', ',
                          plot_cap),
        study = stri_split_fixed(model, 
                                 pattern = '.', 
                                 simplify = TRUE)[, 1],
        progression = stri_split_fixed(model, 
                                       pattern = '.', 
                                       simplify = TRUE)[, 2],
        stage = ifelse(progression == 'early', 'T1 - T2', 'T3 - T4'), 
        x_ax = paste(globals$cohort_lexicon$label[study],
                     stage,
                     sep = ', '))
  
# significance plots -------
  
  insert_msg('Summary bar plots with p values for each gene')
  
  ## summary plot titles
  
  stage_cut$summ_plot_titles <- globals$cxcl_genes %>% 
    map(~paste0('<b><em>', .x, '</em>, survival</b>'))
  
  ## summary plots
  
  stage_cut$summ_plots <-
    list(data = stage_cut$test_results,
         plot_title = stage_cut$summ_plot_titles) %>%
    pmap(plot_signifcant,
         p_variable = 'p_adjusted',
         label_variable = 'x_ax',
         top_significant = 20,
         plot_subtitle = 'Optimal cutoff, Mentel-Henszel test',
         x_lab = expression('-log'[10] * ' pFDR'),
         cust_theme = globals$common_theme) %>% 
    map(~.x + theme(plot.title = element_markdown()))
  
# Kaplan-Meier plots -----
  
  insert_msg('Kaplan-Meier plots')
  
  ## plot titles
  
  stage_cut$plot_titles <- stage_cut$cut_obj %>% 
    map2(., names(.), 
         function(cohort, gene) names(cohort) %>% 
           map(~c(stri_split_fixed(.x, 
                                   pattern = '.', 
                                   simplify = TRUE)[, 1] %>% 
                    exchange(dict = globals$cohort_lexicon, 
                             key = 'cohort', 
                             value = 'label'), 
                 ifelse(stri_detect(.x, fixed = 'early'), 
                        'T1 - T2', 'T3 - T4'))) %>% 
           map(paste, collapse = ' ') %>% 
           map(~paste0('<b><em>', gene, '</em>, ', .x, '</b>')))
  
  
  ## KM plots
  
  stage_cut$km_plots <- 
    list(cut_obj = stage_cut$cut_obj, 
         plot_title = stage_cut$plot_titles) %>% 
    pmap(function(cut_obj, plot_title) list(x = cut_obj, 
                                            title = plot_title) %>% 
           pmap(plot, 
                type = 'km') %>% 
           map(~.x$plot +
                 globals$common_theme + 
                 theme(plot.title = element_markdown())))

  ## appending with the raw and p values
  
  stage_cut$km_plots <- 
    list(x = stage_cut$km_plots, 
         y = stage_cut$test_results) %>%
    pmap(function(x, y) map2(x,
                             y$plot_cap,
                             ~.x + labs(subtitle = .y)))
  
# END --------
  
  insert_tail()