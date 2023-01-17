# Semi-supervised clustering of the tumor samples by lymphokine gene 
# expression.
# The clusters are trained in the TCGA cohort with the PAM/cosine distance 
# algorithm and predicted for the remaining collectives with inverse distance
# weighted kNN classifier

  insert_head()
  
# container ------
  
  gene_clust <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## clustering tables: median centering of the gene expression values
  ## taken over from the clustering development script
  
  gene_clust$genes <- clust_dev$genes
  
  gene_clust$analysis_tbl <-clust_dev$analysis_tbl
  
# semi-supervised clustering -------
  
  insert_msg('Semi-supervised clustering')
  
  ## the training TCGA cohort: the clustering structure was trained before
  ## setting the cluster names
  
  gene_clust$clust_obj$tcga <- clust_dev$algos$pam_cosine
  
  gene_clust$clust_obj$tcga$clust_assignment <- 
    gene_clust$clust_obj$tcga$clust_assignment %>% 
    mutate(clust_id = car::recode(as.character(clust_id), 
                                  "'1' = 'chemox low';  
                               '2' = 'chemox high'"),
           clust_id = factor(clust_id, c('chemox low', 'chemox high')))
  
  ## the test collectives
  
  gene_clust$clust_obj <- 
    gene_clust$analysis_tbl[names(gene_clust$analysis_tbl) != 'tcga'] %>% 
    future_map(~predict.clust_analysis(gene_clust$clust_obj$tcga, 
                                       newdata = .x, 
                                       type = 'propagation', 
                                       kNN = 15, 
                                       simple_vote = FALSE, 
                                       resolve_ties = TRUE), 
               .options = furrr_options(seed = TRUE)) %>% 
    c(gene_clust$clust_obj, .)
  
  ## setting the proper cluster order
  
  for(i in names(gene_clust$clust_obj)) {
    
    gene_clust$clust_obj[[i]]$clust_assignment <- 
      gene_clust$clust_obj[[i]]$clust_assignment %>% 
      mutate(clust_id = factor(clust_id, c('chemox low', 'chemox high')))
    
  }
  
# N numbers --------
  
  insert_msg('N numbers')
  
  gene_clust$n_numbers <- gene_clust$clust_obj %>% 
    map(ngroups)
  
  gene_clust$plot_facets <- gene_clust$n_numbers %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ') %>% 
          set_names(levels(gene_clust$clust_obj[[1]]$clust_assignment$clust_id)))
  
  gene_clust$plot_tags <- gene_clust$plot_facets %>% 
    map(paste, collapse = ', ') %>% 
    map(~paste0('\n', .x))
  
# clustering variances ------
  
  insert_msg('Clustering variances')
  
  gene_clust$variance <- gene_clust$clust_obj %>% 
    map(clustTools::var) %>% 
    map_dbl(~.x$frac_var) %>% 
    compress(names_to = 'cohort', 
             values_to = 'variance') %>% 
    mutate(cohort_lab = exchange(cohort, 
                                 dict = globals$cohort_lexicon, 
                                 key = 'cohort', 
                                 value = 'label'), 
           n = map_dbl(gene_clust$analysis_tbl, nrow), 
           ax_lab = paste(cohort_lab, n, sep = '\nn = '), 
           cohort_type = ifelse(cohort == 'tcga', 
                                'training', 'test'), 
           cohort_type = factor(cohort_type, c('training', 'test')))
  
# Plotting the variances ------
  
  insert_msg('Plotting the variances')
  
  gene_clust$variance_plot <- gene_clust$variance %>% 
    ggplot(aes(x = variance, 
               y = reorder(ax_lab, variance), 
               fill = cohort_type)) +
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_text(aes(label = signif(variance, 2)), 
              color = 'white',
              hjust = 1.4, 
              size = 2.6) + 
    scale_fill_manual(values = c(training = 'coral3', 
                                 test = 'steelblue4'), 
                      labels = c(training = 'PAM/cosine\ntraining cohort', 
                                 test = '15-NN\ntest cohort'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Semi-supervised clustering by lymphokine genes', 
         subtitle = 'Comparison of explained clustering variances', 
         x = 'Explained clustering variance fraction')
  
# Quality control of the clustering structures -------
  
  insert_msg('QC plots')
  
  ## WSS and silhouette curves
  
  gene_clust$training_qc_plots <- gene_clust$clust_obj$tcga %>% 
    plot(type = 'diagnostic', 
         cust_theme = globals$common_theme)
  
  ## distance heat maps
  
  gene_clust$dist_hm <- gene_clust$clust_obj %>% 
    map(plot, 
        type = 'heat_map', 
        cust_theme = globals$common_theme) %>% 
    map2(., globals$cohort_lexicon$label, 
         ~.x + 
          labs(title = .y) + 
           theme(axis.text = element_blank(), 
                 axis.text.x = element_blank(), 
                 axis.title = element_blank(), 
                 axis.line = element_blank(), 
                 axis.ticks = element_blank())) %>% 
    map2(., gene_clust$plot_tags, 
         ~.x + labs(tag = .y))
  
  ## 2D-UMAP
  
  gene_clust$umap_plots <- gene_clust$clust_obj %>% 
    future_map(plot.clust_analysis, 
               type = 'components', 
               cust_theme = globals$common_theme, 
               red_fun = 'umap', 
               with = 'data', 
               kdim = 2, 
               random.state = 123, 
               .options = furrr_options(seed = TRUE, 
                                        packages = c('clustTools'))) %>% 
    map2(., globals$cohort_lexicon$label, 
         ~.x + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           labs(title = .y)) %>% 
    map2(., gene_clust$plot_tags, 
         ~.x + labs(tag = .y))
  
# Heat maps of the clustering features ------
  
  insert_msg('Heat maps of the clustering features')
  
  ## setting the same color scale, squishing to (-3, 3) SDs
  
  gene_clust$feature_hm <- 
    list(x_object = gene_clust$clust_obj, 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle = c('PAM/cosine, training cohort', 
                           rep('15-NN, test cohort', 
                               length(gene_clust$clust_obj) - 1))) %>% 
    pmap(plot_clust_hm, 
         x_lab = 'Cancer sample', 
         cust_theme = globals$common_theme) %>% 
    map2(., gene_clust$plot_facets, 
         ~.x + 
           facet_grid(. ~ sample_clust, 
                      scales = 'free', 
                      space = 'free', 
                      labeller = as_labeller(.y)) + 
           scale_fill_gradient2(low = 'steelblue', 
                                mid = 'black', 
                                high = 'firebrick', 
                                midpoint = 0, 
                                limits = c(-3, 3), 
                                oob = scales::squish, 
                                name = 'Expression Z-score') + 
           scale_y_discrete(limits = rev(gene_clust$genes)) + 
           theme(plot.tag = element_blank(), 
                 axis.text.y = element_text(face = 'italic')))
  
# Tables with normalized expression values of the clustering features -----
  
  insert_msg('Expression tables')
  
  gene_clust$expression_tbl <- gene_clust$clust_obj %>% 
    map(~left_join(set_names(.x$clust_assignment, c('patient_id', 'clust_id')), 
                   rownames_to_column(model.frame(.x), 'patient_id'), 
                   by = 'patient_id'))
  
# Testing for differences in expression between the clusters -----
  
  insert_msg('Testing for differences in gene expression')
  
  ## two-tailed T test
  
  gene_clust$test_results <-  gene_clust$expression_tbl %>% 
    map(~compare_variables(.x, 
                           variables = clust_dev$genes, 
                           split_factor = 'clust_id', 
                           what = 'eff_size', 
                           types = 'cohen_d', 
                           ci = FALSE, 
                           pub_styled = TRUE, 
                           adj_method = 'BH')) %>% 
    map(mutate, 
        plot_cap = paste(eff_size, significance, sep = ', '), 
        plot_lab = paste0('<em>', variable, '</em><br>', significance))
  
  ## plot labels
  
  gene_clust$plot_labs <- gene_clust$test_results %>% 
    map(~set_names(.x$plot_lab, 
                   .x$variable))
  
# Ribbon plots of the clustering features -------
  
  insert_msg('Ribbon plots of the clustering features')
  
  gene_clust$ribbon_plots <- 
    list(data = gene_clust$expression_tbl, 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle = c('PAM/cosine, training cohort', 
                           rep('15-NN, test cohort', 
                               length(gene_clust$clust_obj) - 1))) %>% 
    pmap(draw_stat_panel, 
         variables = gene_clust$genes, 
         split_factor = 'clust_id', 
         stat = 'mean', 
         err_stat = '2se', 
         form = 'line', 
         alpha = 0.25, 
         x_lab = 'mean Z-score, 2 \u00D7 SEM', 
         cust_theme = globals$common_theme) %>% 
    map2(., gene_clust$plot_labs, 
         ~.x + 
           scale_y_discrete(limits = rev(gene_clust$genes), 
                            labels = .y) + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           scale_color_manual(values = globals$cluster_colors, 
                              name = 'Lymphokine cluster') + 
           theme(axis.text.y = element_markdown(), 
                 axis.title.y = element_blank())) %>% 
    map2(., gene_clust$plot_tags, 
         ~.x + labs(tag = .y))

# Fraction of samples assigned to the clusters -------
  
  insert_msg('Fractions of samples assigned to the clusters')
  
  ## plotting data
  
  gene_clust$perc_plot_data <- gene_clust$n_numbers %>% 
    map(mutate, 
        n_total = sum(n), 
        perc = n/n_total * 100) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort_lab = exchange(cohort, 
                                 dict = globals$cohort_lexicon, 
                                 key = 'cohort', 
                                 value = 'label'), 
           cohort_lab = paste(cohort_lab, n_total, sep = '\nn = '))

  ## stack plot
  
  gene_clust$perc_plot <- gene_clust$perc_plot_data %>% 
    ggplot(aes(x = perc, 
               y = cohort_lab, 
               fill = clust_id)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             width = 1.8) + 
    geom_vline(xintercept = 50, 
               linetype = 'dashed') + 
    scale_y_discrete(limits = rev(gene_clust$perc_plot_data$cohort_lab)) + 
    scale_fill_manual(values = globals$cluster_colors, 
                      name = 'Lymphokine cluster') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Distribution of cancer samples', 
         x = '% of cancer samples')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()