# Main figures of the paper

  insert_head()
  
# container -----
  
  figures <- list()
  
# Figure 1: normal vs tumor expression of the lymphokines --------
  
  insert_msg('Figure 1: normal vs tumor expression of the lymphokines')
  
  ## n numbers of the tissue pairs will be provided in the figure legend
  ## effect sizes removed
  
  figures$norm_tumor <- norm_tumor$plots %>% 
    transpose %>% 
    unlist(recursive = FALSE) %>% 
    map(~.x + 
          labs(subtitle = .x$labels$subtitle %>% 
                 stri_replace(regex = ',\\s{1}n.*$', 
                              replacement = '') %>% 
                 stri_replace(regex = '^d.*,\\s{1}', 
                              replacement = '')) + 
          theme(legend.position = 'none')) %>% 
    plot_grid(plotlist = ., 
              ncol = 4, 
              align = 'hv', 
              labels = c('A', '', 'B', '', 
                         'C', '', 'D', '', 
                         'E', '', 'F', '', 
                         'G', '', 'H', ''), 
              label_size = 10) %>% 
    as_figure(label = 'figure_1_normal_tumor', 
              ref_name = 'norm_tumor', 
              caption = paste('Differences in expression of', 
                              'CXCL9/10/11, CXCL13, XCL1 and their', 
                              'cognate receptors between the RCC tumor', 
                              'tissue and normal kidney.'),
              w = 180, 
              h = 220)
  
  
# Figure 2: correlation with infiltration, xcell -------
  
  insert_msg('Figure 2: correlation with infiltration, xcell')
  
  ## force directed plots for the major cohorts
  
  figures$xcell <- xcell_mds$fd_plots[c('tcga', 
                                        'emtab1980', 
                                        'gse73731', 
                                        'gse167093', 
                                        'reca')] %>% 
    map(~.x + 
          scale_alpha_continuous(range = c(0.05, 0.2), 
                                 limits = c(0.1, 0.8), 
                                 oob = scales::squish) + 
          scale_color_manual(labels = c(gene = 'gene', 
                                        cell = 'cell neighbor, \u03C1 > 0.4'), 
                             values = c(gene = 'firebrick4', 
                                        cell = 'steelblue4')))
  
  figures$xcell <- figures$xcell %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(figures$xcell[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2) %>% 
    as_figure(label = 'figure_2_infiltration_xcell', 
              ref_name = 'xcell', 
              caption = paste('Expression of CXCL9/10/11, CXCL13, XCL1', 
                              'and their cognate receptors', 
                              'in RCC samples enriched', 
                              'in CD8^+^ and CD4^+^ T cells,', 
                              'activated myeloid DC, TAM and B cells.'), 
              w = 180, 
              h = 220)
  
# Figure 3: single cell RNAseq results -------
  
  insert_msg('Figure 3: single cell RNAseq results, placeholder')
  
  figures$sc_rnaseq <- ggdraw() %>% 
    as_figure(label = 'figure_3_sc_rnaseq', 
              ref_name = 'sc_rnaseq', 
              caption = paste('Results of scRNAseq analysis, placeholder'), 
              w = 180, 
              h = 180)
  
# Figures 4 and 5: definition and infiltration of the lymphokine clusters ------
  
  insert_msg('Figure 4 and 5: The lymphokine clusters and their infiltration levels')
  
  ## data presented for the major five cohorts
  
  ## Figure 4: definition
  
  figures$lympho_clust <- gene_clust$ribbon_plots[c('tcga', 
                                                    'emtab1980', 
                                                    'gse73731', 
                                                    'gse167093', 
                                                    'reca')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace(.x$labels$tag, 
                                       fixed = '\n', 
                                       replacement = '') %>% 
                 stri_replace_all(fixed = 'chemox ', 
                                  replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(gene_clust$ribbon_plots[[1]])))
  
  ## Figure 5: infiltration by T cells and TAMs
  
  figures$clust_infil <- gene_xcell$plots[c('tcga', 
                                            'emtab1980', 
                                            'gse73731', 
                                            'gse167093', 
                                            'reca')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace_all(.x$labels$subtitle, 
                                           fixed = 'chemox ', 
                                           replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(gene_xcell$plots[[1]])))
  
  ## figure objects
  
  figures[c('lympho_clust', 
            'clust_infil')] <- 
    list(plotlist = figures[c('lympho_clust', 
                              'clust_infil')], 
         ncol = c(3, 2)) %>% 
    pmap(plot_grid, 
         align = 'hv', 
         axis = 'tblr')
  
  figures[c('lympho_clust', 
            'clust_infil')] <- figures[c('lympho_clust', 
                                         'clust_infil')] %>% 
    list(x = ., 
         label = c('figure_4_lymphokine_clusters', 
                   'figure_5_cluster_infiltration'), 
         ref_name = c('lympho_clust', 
                      'clust_infil'), 
         caption = c(paste('Expression of CXCL9/10/11, CXCL13, XCL1', 
                           'and their cognate receptors', 
                           'in the lymphokine clusters.'), 
                     paste('Predicted infiltration of T cells,', 
                           'activated mDC and TAM', 
                           'in the lymphokine clusters.')), 
         h = c(150, 230)) %>% 
    pmap(as_figure, 
         w = 180)

# Figure 6: overall survival in the clusters -------
  
  insert_msg('Figure 6: overall survival in the clusters')
  
  figures$clust_os <- clust_surv$km_plots$os %>% 
    map(~.x + 
          theme(legend.position = 'bottom') + 
          labs(color = NULL)) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    as_figure(label = 'figure_6_cluster_overall_survival', 
              ref_name = 'clust_os', 
              caption = paste('Overall survival in the lymphokine clusters.'), 
              w = 180, 
              h = 220)
  
# Figure 7: biological processes in the clusters -------
  
  insert_msg('Figure 7: biological processes in the clusters')
  
  figures$clust_reactome <- 
    plot_grid(clust_reactome$cmm_volcano_plots$tcga + 
                theme(plot.subtitle = element_blank())) %>% 
    as_figure(label = 'figure_7_cluster_reactome', 
              ref_name = 'clust_reactome', 
              caption = paste('Reactome pathway gene signatures', 
                              'significantly regulated between', 
                              'the lymphokine clusters.'), 
              w = 180, 
              h = 120)
  
# Figure 8: checkpoint genes in the clusters -------
  
  insert_msg('Figure 8: checkpoint genes in the clusters')
  
  figures$checkpoint <- dge_exhaust$checkpoint_plots[c('tcga', 
                                                       'emtab1980', 
                                                       'gse73731', 
                                                       'gse167093', 
                                                       'reca')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace_all(.x$labels$subtitle, 
                                           fixed = 'chemox ', 
                                           replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(dge_exhaust$checkpoint_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_8_clusters_checkpoint', 
              ref_name = 'checkpoint', 
              caption = paste('Differential expression of genes coding', 
                              'for clinically relevant immune checkpoint', 
                              'proteins in the lymphokine clusters.'), 
              w = 180, 
              h = 200)
  
# Figure 9: Metabolism in the clusters ---------
  
  insert_msg('Figure 9: metabolism in the clusters, global')
  
  ## percentages of differentially regulated reactions, TCGA
  ## labeling the pathways of interest with bold
  
  figures$meta_global <- meta_plots$subsystem_plots %>% 
    map(~.x$tcga) %>% 
    map(~.x + 
          scale_y_discrete(labels = metabolism_labeller) + 
          theme(axis.text.y = element_markdown(), 
                plot.title.position = 'plot')) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_9_clusters_metabolism', 
              ref_name = 'meta_global', 
              caption = paste('Differential regulation of metabolic reactions', 
                              'between the lymphocyte clusters.'), 
              w = 180, 
              h = 140)
  
# Figure 10: regulation of the TRP/KYN/QUIN pathway ------
  
  insert_msg('Figure 10: regulation of the TRP/KYN/QUIN pathway genes')
  
  figures$trycats_genes <- meta_genes$trycats_plots[c('tcga', 
                                                      'emtab1980', 
                                                      'gse73731', 
                                                      'gse167093', 
                                                      'reca')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace_all(.x$labels$subtitle, 
                                           fixed = 'chemox ', 
                                           replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(meta_genes$trycats_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(ggdraw() + 
                draw_image('./schemes/trycats.png') + 
                theme(plot.margin = ggplot2::margin(t = 2, l = 10, 
                                                    r = 10, b = 3, 
                                                    unit = 'mm')), 
              ., 
              nrow = 2, 
              rel_heights = c(1, 4.6), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_10_trycats_genes', 
              ref_name = 'trycats_genes', 
              caption = paste('Differential expression of genes involved', 
                              'in the tryptophan/kynurenine/quinolinate', 
                              'metabolism in the lymphokine clusters.'), 
              w = 180, 
              h = 220)
  
# Figure 11: schematic representation of the findings --------
  
  insert_msg('Figure 11: schematic representation of the findings')
  
  figures$summary <- 
    plot_grid(ggdraw() + 
                draw_image('./schemes/lymphokines.png') + 
                theme(plot.margin = ggplot2::margin(t = 3, b = 3, 
                                                    l = 5, r = 5, 
                                                    unit = 'mm')), 
              ggdraw() + 
                draw_image('./schemes/immunosuppression.png') + 
                theme(plot.margin = ggplot2::margin(t = 3, b = 3, 
                                                    l = 5, r = 5, 
                                                    unit = 'mm')), 
              ncol = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_11_summary', 
              ref_name = 'summary', 
              caption = paste('Recruitment of CD8^+^ T cells via CXCL9/10/11,', 
                              'CXCL13 and XCL1 and suppression of CD8^+^', 
                              'T cells in the microenvironment of RCC.'), 
              w = 180, 
              h = 80)

# saving the figures ------
  
  insert_msg('Saving the figures')
  
  figures %>% 
    walk(pickle, 
         path = './paper/figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
  
# END ------
  
  insert_tail()