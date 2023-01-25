# Comparison of the xCell infiltration levels between the chemokine clusters
#
# Done with Wilcoxon test with r effect size statistic

  insert_head()
  
# container -----
  
  gene_xcell <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  gene_xcell$variables <- globals$xcell_lexicon$variable
  
  ## analysis tables
  
  gene_xcell$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('gene_xcell$analysis_tbl$', i, 
                           ' <- ', i, '$xcell')))

  }
  
  gene_xcell$analysis_tbl <- gene_clust$clust_obj %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    map2(gene_xcell$analysis_tbl, 
         left_join, by = 'patient_id')
  
# Descriptive stats -------
  
  insert_msg('Descriptive stats')
  
  gene_xcell$desc_stats <- gene_xcell$analysis_tbl %>% 
    future_map(explore, 
               variables = gene_xcell$variables, 
               split_factor = 'clust_id', 
               what = 'table', 
               pub_styled = TRUE, 
               .options = furrr_options(seed = TRUE)) %>% 
    map(reduce, 
        left_join, 
        by = 'variable') %>% 
    map(set_names, 
        c('variable', levels(gene_xcell$analysis_tbl[[1]]$clust_id)))
  
# testing -------
  
  insert_msg('Testing')
  
  gene_xcell$test_results <- gene_xcell$analysis_tbl %>% 
    future_map(compare_variables, 
               variables = gene_xcell$variables, 
               split_factor = 'clust_id', 
               what = 'eff_size', 
               types = 'wilcoxon_r', 
               pub_styled = TRUE, 
               ci = FALSE, 
               exact = FALSE, 
               adj_method = 'BH', 
               .options = furrr_options(seed = TRUE)) %>% 
    map(mutate, 
        var_lab = exchange(variable, 
                           dict = globals$xcell_lexicon, 
                           key = 'variable', 
                           value = 'label'), 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        ax_lab = paste(var_lab, significance, sep = '<br>'), 
        ax_lab = ifelse(p_adjusted < 0.05, 
                        paste0('<b>', ax_lab, '</b>'), 
                        ax_lab))
  
# cell populations differing between the clusters in all cohorts ------
  
  insert_msg('Significant and near significant cell types')
  
  ## defined as p raw < 0.05 in all cohorts except the smallest CM10 cohort
  
  gene_xcell$significant <- 
    gene_xcell$test_results[names(gene_xcell$test_results) != 'cm10'] %>% 
    map(filter, p_value < 0.05) %>% 
    map(~.x$variable) %>% 
    reduce(intersect)
  
# Plotting the selected cell type levels in box plots -------
  
  insert_msg('Box plots')
  
  ## plots for CD8+ central memory T cells, CD8+ effector memory T cells
  ## CD4+ memory, CD+ effector memory, activated mDC, M1 and M2 TAMs
  
  ## cell types
  
  gene_xcell$selected_cells <- 
    c('T cell CD8+ central memory', 
      'T cell CD8+ effector memory', 
      'T cell CD4+ memory', 
      'T cell CD4+ effector memory', 
      'B cell', 
      'Myeloid dendritic cell activated', 
      'Macrophage M1', 
      'Macrophage M2')
  
  ## axis labels
  
  gene_xcell$axis_labs <- gene_xcell$test_results %>% 
    map(~set_names(.x$ax_lab, .x$variable))
  
  ## plotting tables with normalized infiltration levels
  
  gene_xcell$plot_tbl <- gene_xcell$analysis_tbl %>% 
    map(select, 
        patient_id, 
        clust_id, 
        all_of(gene_xcell$selected_cells)) %>%
    map(~map_dfc(.x, function(x) if(is.numeric(x)) scale(x, center = median(x))[, 1] else x))
  
  ## n numbers to be presented in the plot captions
  
  gene_xcell$plot_caps <- gene_xcell$plot_tbl %>% 
    map(count, clust_id) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ')) %>% 
    map(paste, collapse = ', ')
  
  ## plotting
  
  gene_xcell$plots <- list(data = gene_xcell$plot_tbl, 
                           plot_title = globals$cohort_lexicon$label, 
                           plot_subtitle = gene_xcell$plot_caps) %>% 
    pmap(draw_violin_panel, 
         variables = gene_xcell$selected_cells, 
         split_factor = 'clust_id', 
         distr_geom = 'box', 
         point_hjitter = 0, 
         point_wjitter = 0.1, 
         point_alpha = 0.1, 
         point_size = 1.5, 
         cust_theme = globals$common_theme, 
         x_lab = 'Z score, xCell', 
         outlier.color = NA) %>% 
    map2(., gene_xcell$axis_labs, 
         ~.x + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           scale_y_discrete(labels = .y, 
                            limits = rev(gene_xcell$selected_cells)) + 
           theme(axis.text.y = element_markdown()))
  
# END ------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()
