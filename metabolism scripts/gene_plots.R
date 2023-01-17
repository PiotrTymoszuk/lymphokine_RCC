# Plots for selected metabolic pathway genes: 
# OxPhos, TCA, FAO, glycolysis/glyconeogenesis, tryptophan metabolism 
# and IDO/TRACATS

  insert_head()
  
# container ------
  
  meta_genes <- list()
  
# Volcano and top gene plots for FAO reactions -----
  
  insert_msg('FAO reaction genes')
  
  meta_genes$fao_genes <- sig$recon_lexicon %>% 
    filter(variable == 'fatty_acid_oxidation') %>% 
    .$gene_symbol %>% 
    reduce(union)
  
  meta_genes$fao_genes <- meta_genes$fao_genes[!is.na(meta_genes$fao_genes)]
  
  meta_genes$fao_volcano <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$fao_genes), 
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
  
  meta_genes$fao_top <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$fao_genes), 
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
  
# TCA genes -------
  
  insert_msg('TCA reaction genes')
  
  meta_genes$tca_genes <- sig$recon_lexicon %>% 
    filter(variable == 'citric_acid_cycle') %>% 
    .$gene_symbol %>% 
    reduce(union)
  
  meta_genes$tca_genes <- meta_genes$tca_genes[!is.na(meta_genes$tca_genes)]
  
  meta_genes$tca_volcano <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$tca_genes), 
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
  
  meta_genes$tca_top <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$tca_genes), 
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
  
# OxPhos genes ------
  
  insert_msg('OxPhos reaction genes')

  meta_genes$oxphos_genes <- sig$recon_lexicon %>% 
    filter(variable == 'oxidative_phosphorylation') %>% 
    .$gene_symbol %>% 
    reduce(union)
  
  meta_genes$oxphos_genes <- 
    meta_genes$oxphos_genes[!is.na(meta_genes$oxphos_genes)]
  
  meta_genes$oxphos_volcano <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$oxphos_genes), 
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
  
  meta_genes$oxphos_top <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$oxphos_genes), 
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
  

# TRP metabolism genes ----
  
  insert_msg('TRP metabolism genes')
  
  meta_genes$trp_genes <- sig$recon_lexicon %>% 
    filter(variable == 'tryptophan_metabolism') %>% 
    .$gene_symbol %>% 
    reduce(union)
  
  meta_genes$trp_genes <- 
    meta_genes$trp_genes[!is.na(meta_genes$trp_genes)]
  
  meta_genes$trp_volcano <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$trp_genes), 
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
  
  meta_genes$trp_top <- 
    list(data = clust_dge$test_results %>% 
           map(filter, gene_symbol %in% meta_genes$trp_genes), 
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
  
# TRYCATS genes -------
  
  insert_msg('Trycats genes')
  
  ## genes: BH4 biosynthesis and TRP decay to KYN/QUIN
  
  meta_genes$trycats_genes <- 
    c('IDO1', 'IDO2', 'TDO2',
      'KYAT1', 'KYNU', 'KMO', 'HAAO')
  
  ## testing results and plot labels
  
  meta_genes$trycats_test_results <- clust_dge$test_results %>% 
    map(filter, gene_symbol %in% meta_genes$trycats_genes) %>% 
    map(mutate, 
        significance = ifelse(p_adjusted >= 0.05, 
                              paste0('ns (p = ', signif(p_adjusted, 2), ')'), 
                              paste('p =', signif(p_adjusted, 2))), 
        ax_lab = paste0('<em>', gene_symbol, '</em><br>', significance), 
        ax_lab = ifelse(p_adjusted < 0.05, 
                        paste0('<b>', ax_lab, '</b>'), 
                        ax_lab))
  
  meta_genes$trycats_ax_labs <- meta_genes$trycats_test_results %>% 
    map(~set_names(.x$ax_lab, 
                   .x$gene_symbol))
  
  ## tables with normalized expression values
  
  meta_genes$trycats_expression <- clust_dge$analysis_tbl %>% 
    map(select, clust_id, any_of(meta_genes$trycats_genes)) %>% 
    map(~map_dfc(.x, function(x) if(is.numeric(x)) scale(x)[, 1] else x))
  
  ## plotting variables
  
  meta_genes$trycats_variables <- meta_genes$trycats_expression %>% 
    map(names) %>% 
    map(~.x[.x %in% meta_genes$trycats_genes])
  
  ## plot captions with n numbers of the clusters
  
  meta_genes$trycats_captions <- meta_genes$trycats_expression %>% 
    map(count, clust_id) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ')) %>% 
    map(paste, collapse = ', ')
  
  ## panels of box plots
  
  meta_genes$trycats_plots <- 
    list(data = meta_genes$trycats_expression, 
         variables = meta_genes$trycats_variables, 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle = meta_genes$trycats_captions) %>% 
    pmap(draw_violin_panel, 
         split_factor = 'clust_id', 
         distr_geom = 'box', 
         point_alpha = 0.1, 
         point_hjitter = 0, 
         point_wjitter = 0.15,
         point_size = 1.5, 
         x_lab = 'expression Z score', 
         cust_theme = globals$common_theme, 
         outlier.color = NA)
  
  meta_genes$trycats_plots <- 
    list(x = meta_genes$trycats_plots, 
         y = meta_genes$trycats_variables, 
         z = meta_genes$trycats_ax_labs) %>% 
    pmap(function(x, y, z) x + 
           scale_y_discrete(limits = rev(y), 
                            labels = z) + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           theme(axis.title.y = element_blank(), 
                 axis.text.y = element_markdown()))

# END -------
  
  insert_tail()