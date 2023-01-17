# Regulation of dge_exhaustion genes

  insert_head()
  
# container ----
  
  dge_exhaust <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## signatures of dge_exhaustion
  ## sources:
  ## https://doi.org/10.1016/j.ebiom.2022.104207
  ## https://doi.org/10.1158/1078-0432.CCR-22-0275
  ## https://doi.org/10.1158/1078-0432.CCR-17-1846
  ## https://doi.org/10.1002/path.5435
  ## There is a discrepant gene VSIR (Entrez ID: 64115)
  ## which has not official symbol in the org.Hs database
  ## but corresponds to the C10orf54 item
  
  dge_exhaust$genes <- c('BTLA',
                         'CTLA4',
                         'CD274',
                         'PDCD1',
                         'PDCD1LG2',
                         'HAVCR2',
                         'IDO1',
                         'SIGLEC7',
                         'TIGIT',
                         'LAG3',
                         'VSIR',
                         'CXCR4',
                         'CD160',
                         'ENTPD1',
                         'CD244',
                         'EOMES',
                         'CD274',
                         'CD200',
                         'TOX',
                         'TOX2',
                         'C10orf54')
  
  dge_exhaust$genes <-
    read_excel('./data/literature signatures/path5435-sup-0002-supptables.xlsx',
               sheet = 'S5',
               skip = 3) %>%
    .$gene %>%
    c(dge_exhaust$gene) %>%
    unique
  
  dge_exhaust$genes <- dge_exhaust$genes[dge_exhaust$genes != 'VSIR']
  
  ## gene labels
  
  dge_exhaust$genes <- ifelse(dge_exhaust$genes == 'C10orf54',
                              'VSIR', dge_exhaust$genes) %>%
    set_names(dge_exhaust$genes) %>%
    compress(names_to = 'gene_symbol',
             values_to = 'gene_label') %>%
    filter(!duplicated(gene_label))
  
# Analysis tables -------
  
  insert_msg('Analysis tables')
  
  dge_exhaust$analysis_tbl <- clust_dge$test_results %>% 
    map(filter, 
        gene_symbol %in% dge_exhaust$genes$gene_symbol)
  
  ## common regulated genes (all cohorts except of the smallest CM10)
  
  dge_exhaust$common_up <- 
    dge_exhaust$analysis_tbl[names(dge_exhaust$analysis_tbl) != 'cm10'] %>% 
    map(filter, regulation == 'upregulated') %>% 
    map(~.x$gene_symbol) %>%
    reduce(intersect)
  
# Volcano plots -------
  
  insert_msg('Volcano plots')
  
  dge_exhaust$volcano_plots <- 
    list(data = dge_exhaust$analysis_tbl, 
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
  
# Top regulation plots ------
  
  insert_msg('Top regulation plots')
  
  dge_exhaust$top_plots <- 
    list(data = dge_exhaust$analysis_tbl %>% 
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
                               'and downregulated genes, T cell exhaustion'), 
         x_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         cust_theme = globals$common_theme) %>% 
    map(~.x + 
          theme(axis.text.y = element_text(face = 'italic')))
  
# Plots for the clinically relevant checkpoint genes ------
  
  insert_msg('Clinically relevant genes')
  
  ## genes
  
  dge_exhaust$checkpoint_genes <- 
    c('CD274', 'PDCD1LG2', 'PDCD1', 'CTLA4', 'HAVCR2')
  
  ## normalized expression data
  
  dge_exhaust$expression <- clust_dge$analysis_tbl %>% 
    map(select, 
        clust_id, 
        any_of(dge_exhaust$checkpoint_genes)) %>% 
    map(~map_dfc(.x, function(x) if(is.numeric(x)) scale(x)[, 1] else x))
  
  ## plot captions with the n numbers
  
  dge_exhaust$plot_captions <- dge_exhaust$expression %>% 
    map(count, clust_id) %>% 
    map(~map2_chr(.x[[1]], .x[[2]], 
                  paste, sep = ': n = ')) %>% 
    map(paste, collapse = ', ')
  
  ## gene labs with significance of regulation between the clusters

  dge_exhaust$x_axes <- dge_exhaust$analysis_tbl %>% 
    map(filter, 
        gene_symbol %in% dge_exhaust$checkpoint_genes) %>% 
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
  
  dge_exhaust$checkpoint_plots <- 
    list(data = dge_exhaust$expression, 
         variables = dge_exhaust$expression %>% 
           map(names)%>% 
           map(~.x[.x %in% dge_exhaust$checkpoint_genes]), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle = dge_exhaust$plot_captions) %>% 
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
  
  dge_exhaust$checkpoint_plots <- 
    list(x = dge_exhaust$checkpoint_plots, 
         y = dge_exhaust$expression %>% 
           map(names)%>% 
           map(~.x[.x %in% dge_exhaust$checkpoint_genes]), 
         z = dge_exhaust$x_axes) %>% 
    pmap(function(x, y, z) x + 
           scale_y_discrete(limits = y, 
                            labels = z) + 
           scale_fill_manual(values = globals$cluster_colors, 
                             name = 'Lymphokine cluster') + 
           theme(axis.text.y = element_markdown()))
# END -----
  
  insert_tail()