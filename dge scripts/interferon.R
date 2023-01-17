# Differential expression of interferon signaling genes

  insert_head()
  
# container -------
  
  dge_ifn <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## genes of interest: extracted from the IFN gamma, alpha and beta
  ## Reactome signatures
  
  dge_ifn$genes <- 
    c('./data/literature signatures/Participating Molecules [R-HSA-877300].tsv', 
      './data/literature signatures/Participating Molecules [R-HSA-909733].tsv') %>% 
    map(read_tsv) %>%  
    map(mutate, 
        gene_symbol = stri_split_fixed(MoleculeName, 
                                       pattern = ' ', 
                                       simplify = TRUE)[, 2]) %>% 
    map(~.x$gene_symbol) %>% 
    map(unique) %>% 
    reduce(union)
  
# Analysis tables -------
  
  insert_msg('Analysis tables')
  
  dge_ifn$analysis_tbl <- clust_dge$test_results %>% 
    map(filter, 
        gene_symbol %in% dge_ifn$genes)
  
# numbers of significantly up- and downregulated genes -------
  
  insert_msg('Numbers of significantly regulated genes')

  dge_ifn$regulation_counts <- dge_ifn$analysis_tbl %>% 
    map(filter, regulation != 'ns') %>% 
    map(count, regulation, .drop = FALSE) %>% 
    compress(names_to = 'cohort') %>% 
    filter(regulation != 'ns')
  
# common regulated genes (all cohorts except of the smallest CM10) -----
  
  insert_msg('Common regulated genes')
  
  dge_ifn$common_up <- 
    dge_ifn$analysis_tbl[names(dge_ifn$analysis_tbl) != 'cm10'] %>% 
    map(filter, regulation == 'upregulated') %>% 
    map(~.x$gene_symbol) %>%
    reduce(intersect)
  
  dge_ifn$common_down <- 
    dge_ifn$analysis_tbl[names(dge_ifn$analysis_tbl) != 'cm10'] %>% 
    map(filter, regulation == 'downregulated') %>% 
    map(~.x$gene_symbol) %>%
    reduce(intersect)
  
# Volcano plots -------
  
  insert_msg('Volcano plots')
  
  dge_ifn$volcano_plots <- 
    list(data = dge_ifn$analysis_tbl, 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_volcano, 
         regulation_variable = 'estimate', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = log2(1.25), 
         x_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         y_lab = expression('-log'[10] * ' pFDR'), 
         top_significant = 30, 
         label_variable = 'gene_symbol', 
         label_type = 'text', 
         txt_size = 2.4, 
         txt_face = 3, 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          labs(subtitle = .x$labels$tag) + 
          theme(plot.tag = element_blank()))
  
# Top regulation plots ------
  
  insert_msg('Top regulation plots')
  
  dge_ifn$top_plots <- 
    list(data = dge_ifn$analysis_tbl %>% 
           map(filter, regulation != 'ns'), 
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
                               'and downregulated genes, IFN signaling'), 
         x_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          theme(axis.text.y = element_text(face = 'italic')))
  
# Plots for the common regulated pathway genes ------
  
  insert_msg('Common regulated genes')

  ## normalized expression data
  
  dge_ifn$expression <- clust_dge$analysis_tbl %>% 
    map(select, 
        clust_id, 
        any_of(dge_ifn$common_up)) %>% 
    map(~map_dfc(.x, function(x) if(is.numeric(x)) scale(x)[, 1] else x))
  
  ## plot captions with the n numbers
  
  dge_ifn$plot_captions <- dge_ifn$expression %>% 
    map(count, clust_id) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ')) %>% 
    map(paste, collapse = ', ')
  
  ## gene labs with significance of regulation between the clusters
  
  dge_ifn$x_axes <- dge_ifn$analysis_tbl %>% 
    map(filter, 
        gene_symbol %in% dge_ifn$common_up) %>% 
    map(mutate, 
        significance = ifelse(p_adjusted < 0.05, 
                              paste('p =', signif(p_adjusted, 2)), 
                              paste0('ns (p = ', signif(p_adjusted, 2), ')')), 
        x_axes = paste0('<em>', gene_symbol, '</em><br>', significance), 
        x_axes = ifelse(p_adjusted < 0.05, 
                        paste0('<b>', x_axes, '</b>'), 
                        x_axes)) %>% 
    map(~set_names(.x$x_axes, .x$gene_symbol))
  
  ## plots 
  
  dge_ifn$common_up_plots <- 
    list(data = dge_ifn$expression, 
         variables = dge_ifn$expression %>% 
           map(names)%>% 
           map(~.x[.x %in% dge_ifn$common_up]), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle = dge_ifn$plot_captions) %>% 
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
  
  dge_ifn$common_up_plots <- 
    list(x = dge_ifn$common_up_plots, 
         y = dge_ifn$expression %>% 
           map(names)%>% 
           map(~.x[.x %in% dge_ifn$common_up]), 
         z = dge_ifn$x_axes) %>% 
    pmap(function(x, y, z) x + 
           scale_y_discrete(limits = y, 
                            labels = z) + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           theme(axis.text.y = element_markdown()))
  
# END -----
  
  insert_tail()