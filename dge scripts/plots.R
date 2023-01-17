# Volcano plots and top regulated gene plots for the differentially 
# regulated genes in the lymphokine clusters.

  insert_head()
  
# container -------
  
  dge_plots <- list()
  
# gene number plots -------
  
  insert_msg('Gene number plots')
  
  ## regulated gene numbers
  
  dge_plots$gene_numbers <- clust_dge$significant %>% 
    map(count, regulation) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(n_recoded = ifelse(regulation == 'downregulated', 
                              -n, n))
  
  ## gene number plot
  
  dge_plots$gene_number_plot <- dge_plots$gene_numbers %>% 
    ggplot(aes(x = n_recoded, 
               y = reorder(cohort, n),
               fill = regulation)) + 
    geom_bar(stat = 'identity', 
             color = 'black') + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    geom_text(aes(label = n, 
                  color = regulation, 
                  x = ifelse(n_recoded < 0, 
                             n_recoded - 100, 
                             n_recoded + 100)), 
              size = 2.6, 
              hjust = 0.5) + 
    scale_y_discrete(labels = globals$cohort_lexicon$label) + 
    scale_fill_manual(values = c(upregulated = 'firebrick4', 
                                 downregulated = 'steelblue4'), 
                      name = 'Regulation\nchemox high vs low') + 
    scale_color_manual(values = c(upregulated = 'firebrick4', 
                                  downregulated = 'steelblue4'), 
                       name = 'Regulation\nchemox high vs low') + 
    scale_x_continuous(labels = function(x) abs(x), 
                       limits = c(-750, 2800), 
                       breaks = seq(-500, 2500, by = 500)) + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Differential gene expression', 
         x = '# regulated genes, chemox high vs low')
  
# Volcano plots ------
  
  insert_msg('Volcano plots')
  
  dge_plots$volcano_plots <- 
    list(data = clust_dge$test_results, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = log2(1.25), 
         x_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         y_lab = expression('-log'[10] * ' pFDR'), 
         top_significant = 20, 
         label_variable = 'gene_symbol', 
         label_type = 'text', 
         txt_size = 2.4, 
         txt_face = 3, 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          labs(subtitle = .x$labels$tag) + 
          theme(plot.tag = element_blank()))
  
# Top regulated plots -------
  
  insert_msg('Top 20 regulated plots')
  
  dge_plots$top_plots <- 
    list(data = clust_dge$significant, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_top, 
         regulation_variable = 'estimate', 
         label_variable = 'gene_symbol', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         lower_ci_variable = 'lower_ci', 
         upper_ci_variable = 'upper_ci', 
         top_regulated = 20, 
         plot_subtitle = paste('Top 20 stongest differentially up-', 
                               'and downregulated genes'), 
         x_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          theme(axis.text.y = element_text(face = 'italic')))
  
# END --------
  
  insert_tail()