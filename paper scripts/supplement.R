# Supplementary figures 

  insert_head()
  
# container ------
  
  suppl_figures <- list()
  
# Figure S1: chemokine screening -------
  
  insert_msg('Figure S1: chemokine screening')
  
  suppl_figures$chemo_screen <- 
    plot_grid(plot_grid(screen$norm_tumor_volcano + 
                          theme(legend.position = 'none', 
                                plot.tag = element_blank()), 
                        plot_grid(get_legend(screen$norm_tumor_volcano), 
                                  ggdraw() + 
                                    draw_text(screen$norm_tumor_volcano$labels$tag, 
                                              size = 8, 
                                              hjust = 0, 
                                              x = 0.3), 
                                  nrow = 2), 
                        ncol = 2, 
                        rel_widths = c(0.7, 0.3)), 
              screen$corr_forest, 
              nrow = 2,
              rel_heights = c(0.4, 0.6), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s1_chemokine_screening', 
              ref = 'chemo_screen', 
              caption = paste('Screening for tumor-expressed chemokines', 
                              'involved in CD8^+^ T cell recruitment in', 
                              'the TCGA KIRC cohort.'), 
              w = 180, 
              h = 210)
  
# Figure S2: correlation with infiltration, networks, xcell -------
  
  insert_msg('Figure S2: correlation with infiltration, xcell')
  
  ## force directed plots for the major cohorts
  
  suppl_figures$xcell_net <- xcell_mds$fd_plots[c('cm10', 
                                                  'cm25ev', 
                                                  'cm25ni')] %>% 
    map(~.x + 
          scale_alpha_continuous(range = c(0.1, 0.25), 
                                 limits = c(0.1, 0.8), 
                                 oob = scales::squish) + 
          scale_color_manual(labels = c(gene = 'gene', 
                                        cell = 'cell neighbor, \u03C1 > 0.4'), 
                             values = c(gene = 'firebrick4', 
                                        cell = 'steelblue4')))
  
  suppl_figures$xcell_net <- suppl_figures$xcell_net %>% 
    map(~.x + theme(legend.position = 'none')) %>% 
    c(list(get_legend(suppl_figures$xcell_net[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2) %>% 
    as_figure(label = 'figure_s2_infiltration_xcell', 
              ref_name = 'xcell_net', 
              caption = paste('Expression of CXCL9/10/11, CXCL13, XCL1', 
                              'and their cognate receptors', 
                              'in RCC samples enriched', 
                              'in CD8^+^ and CD4^+^ T cells,', 
                              'activated myeloid DC, TAM and B cells,', 
                              'CheckMate cohorts.'), 
              w = 180, 
              h = 150)
  
# Figures S3 - S4: infiltration, xCell --------
  
  insert_msg('Figure S3 - S4: infiltration xCell')
  
  suppl_figures[c('xcell1', 'xcell2')] <- 
    list(xcell$heat_bubble[c('tcga', 'emtab1980', 'gse73731', 'gse167093')], 
         xcell$heat_bubble[c('reca', 'cm10', 'cm25ev', 'cm25ni')]) %>% 
    map(~map(.x, ~.x + 
               theme(legend.position = 'none', 
                     axis.text = element_text(size = 6), 
                     plot.margin = ggplot2::margin(t = 2, r = 1, 
                                                   b = 2, l = 1, 
                                                   unit = 'mm')))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = 'hv') %>% 
          plot_grid(get_legend(xcell$heat_bubble[[1]] + 
                                 theme(legend.position = 'bottom', 
                                       legend.text = element_text(size = 6))), 
                    nrow = 2, 
                    rel_heights = c(0.95, 0.05)))
  
  suppl_figures[c('xcell1', 'xcell2')] <- 
    suppl_figures[c('xcell1', 'xcell2')] %>% 
    list(x = ., 
         label = c('figure_s3_xcell', 
                   'figure_s4_xcell'), 
         ref_name = c('xcell1', 'xcell2'), 
         caption = c(paste('Correlation of expression of CXCL9/10/11,', 
                           'CXCL13, XCL1 and their cognate receptors with', 
                           'non-malingant cell population content predicted', 
                           'by the xCell algorithm in the', 
                           'TCGA KIRC, E-MTAB 1980, GSE73731', 
                           'and GSE167093 cohorts.'), 
                     paste('Correlation of expression of CXCL9/10/11,', 
                           'CXCL13, XCL1 and their cognate receptors with', 
                           'non-malingant cell population content predicted', 
                           'by the xCell algorithm in the', 
                           'RECA-EU, CheckMate010 and', 
                           'CheckMate025 cohorts.'))) %>% 
    pmap(as_figure, 
         w = 180, 
         h = 230)
  
# Figures S5 - S6: infiltration, quanTIseq -------
  
  insert_msg('Figures S5 - S6: infiltration, quanTIseq')
  
  suppl_figures[c('quantiseq1', 'quantiseq2')] <- 
    list(quantiseq$heat_bubble[c('tcga', 'emtab1980', 'gse73731', 'gse167093')], 
         quantiseq$heat_bubble[c('reca', 'cm10', 'cm25ev', 'cm25ni')]) %>% 
    map(~map(.x, 
             ~.x + 
               theme(legend.position = 'none', 
                     axis.text = element_text(size = 6), 
                     plot.margin = ggplot2::margin(t = 2, r = 1, 
                                                   b = 2, l = 1, 
                                                   unit = 'mm')))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = 'hv') %>% 
          plot_grid(get_legend(xcell$heat_bubble[[1]] + 
                                 theme(legend.position = 'bottom', 
                                       legend.text = element_text(size = 6))), 
                    nrow = 2, 
                    rel_heights = c(0.95, 0.05)))
  
  suppl_figures[c('quantiseq1', 'quantiseq2')] <- 
    suppl_figures[c('quantiseq1', 'quantiseq2')] %>% 
    list(x = ., 
         label = c('figure_s5_quantiseq', 
                   'figure_s6_quantiseq'), 
         ref_name = c('quantiseq1', 'quantiseq2'), 
         caption = c(paste('Correlation of expression of CXCL9/10/11,', 
                           'CXCL13, XCL1 and their cognate receptors with', 
                           'immune cell population content predicted', 
                           'by the QuanTIses algorithm in the', 
                           'TCGA KIRC, E-MTAB 1980, GSE73731', 
                           'and GSE167093 cohorts.'), 
                     paste('Correlation of expression of CXCL9/10/11,', 
                           'CXCL13, XCL1 and their cognate receptors with', 
                           'immune cell population content predicted', 
                           'by the QuanTIseq algorithm in the', 
                           'RECA-EU, CheckMate010 and', 
                           'CheckMate025 cohorts.'))) %>% 
    pmap(as_figure, 
         w = 180, 
         h = 200)
  
# Figure S7: single cell RNA seq -------
  
  insert_msg('Figure S7: scRNAseq')
  
  suppl_figures$scrna_seq <- 
    plot_grid(ggdraw() + 
                draw_image('./schemes/sc_rna_seq.png')) %>% 
    as_figure(label = 'figure_s7_sc_rna_seq', 
              ref_name = 'scrna_seq', 
              caption = 'Method and sample flowchart of scRNA dataset.', 
              w = 180, 
              h = 4800/5417 * 180)
  
# Figure S8: cluster training -------
  
  insert_msg('Figure S8: training of the clusters')
  
  suppl_figures$clust_training <- 
    plot_grid(plot_grid(clust_dev$result_plot + 
                          theme(legend.position = 'none'), 
                        plot_grid(get_legend(clust_dev$result_plot), 
                                  ggdraw() + 
                                    draw_text(gene_clust$training_qc_plots$wss$labels$tag, 
                                              size = 8, 
                                              hjust = 0, 
                                              x = 0.1), 
                                  nrow = 2), 
                        ncol = 2, 
                        rel_widths = c(0.8, 0.2)), 
              plot_grid(plotlist = gene_clust$training_qc_plots %>% 
                          map(~.x + 
                                labs(title = 'Cluster number, TCGA KIRC') + 
                                theme(plot.tag = element_blank())), 
                        ncol = 2, 
                        align = 'hv'), 
              plot_grid(gene_clust$dist_hm$tcga + 
                          theme(legend.position = 'none', 
                                plot.tag = element_blank()),
                        plot_grid(get_legend(gene_clust$dist_hm$tcga + 
                                               theme(legend.position = 'bottom')), 
                                  ggdraw() + 
                                    draw_text(gene_clust$dist_hm$tcga$labels$tag %>% 
                                                stri_replace_all(fixed = ', ', 
                                                                 replacement = '\n'), 
                                              size = 8, 
                                              hjust = 0, 
                                              x = 0.3), 
                                  nrow = 2), 
                        ncol = 2), 
              nrow = 3, 
              rel_heights = c(1.2, 1, 1), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s8_lymphokine_cluster_training', 
              ref_name = 'clust_training', 
              caption = paste('Definition of the lymphokine', 
                              'clusters in the TCGA KIRC cohort.'), 
              w = 180, 
              h = 230)
  
# Figure S9: semi-supervised clustering -------
  
  insert_msg('Figure S9: semi-supervised clustering')
  
  ## upper panel with clustering variances
  ## and fractions of cancer samples in the clusters
  
  suppl_figures$semi_clust$upper_panel <- 
    plot_grid(gene_clust$variance_plot + 
                labs('Explanatory performance') + 
                theme(plot.subtitle = element_blank(), 
                      legend.position = 'bottom', 
                      legend.title = element_blank()), 
              gene_clust$perc_plot + 
                theme(legend.position = 'bottom', 
                      legend.title = element_blank()), 
              ncol = 2, 
              align = 'hv')

  ## lower panel with the ribbon plots of gene expression
  ## for the CheckMate cohorts
  
  suppl_figures$semi_clust$bottom_panel <- 
    gene_clust$ribbon_plots[c('cm10', 'cm25ev', 'cm25ni')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace(.x$labels$tag, 
                                       fixed = '\n', 
                                       replacement = '') %>% 
                 stri_replace_all(fixed = 'chemox ', 
                                  replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                legend.title = element_blank(), 
                plot.margin = ggplot2::margin(t = 4, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(get_legend(gene_clust$ribbon_plots[[1]] + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.9, 0.1))
  
  ## the entire figure
  
  suppl_figures$semi_clust <- 
    plot_grid(suppl_figures$semi_clust$upper_panel, 
              suppl_figures$semi_clust$bottom_panel, 
              nrow = 2, 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s9_semisupervised_clustering', 
              ref_name = 'semi_clust', 
              caption = paste('Semi-supervised clustering of cancer', 
                              'samples by the lymphokine gene expression.', 
                              'Expression of CXCL9/10/11, CXCL13, XCL1', 
                              'and their cognate receptors', 
                              'in the lymphokine clusters, checkMate cohorts.'), 
              w = 180, 
              h = 180)
  
# Figure S10: infiltration in the checkmate cohorts ------
  
  insert_msg('Figure S10: infiltration in the CheckMate cohorts')
  
  suppl_figures$clust_infil <- gene_xcell$plots[c('cm10', 
                                                  'cm25ev', 
                                                  'cm25ni')] %>% 
    map(~.x + 
          labs(subtitle = stri_replace_all(.x$labels$subtitle, 
                                           fixed = 'chemox ', 
                                           replacement = '')) + 
          theme(legend.position = 'none', 
                plot.tag = element_blank(), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(gene_xcell$plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s10_cluster_infiltration', 
              ref_name = 'clust_infil', 
              caption = paste('Predicted infiltration of T cells,', 
                              'activated mDC and TAM', 
                              'in the lymphokine clusters, CheckMate cohorts.'), 
              w = 180, 
              h = 155)
  
# Figure S11: tumor stage and MSKCC risk group in the lymphokine clusters -----
  
  insert_msg('Figure S11: tumor stage and risk in the clusters')
  
  ## upper panel: staging
  
  suppl_figures$stage$upper_panel <- 
    clust_clinic$plots[c('tcga', 
                         'emtab1980', 
                         'gse73731', 
                         'gse167093')] %>% 
    map2(., globals$cohort_lexicon$label[c('tcga', 
                                           'emtab1980', 
                                           'gse73731', 
                                           'gse167093')], 
         ~.x$pt_stage + 
          theme(legend.position = 'none', 
                axis.title.x = element_blank()) + 
          labs(title = .y, 
               subtitle = .x$pt_stage$labels$subtitle %>% 
                 stri_replace(regex = '^V.*,\\s{1}', 
                              replacement = ''))) %>% 
    c(list(get_legend(clust_clinic$plots$tcga$pt_stage))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv')
  
  ## bottom panel: MSKCC
  
  suppl_figures$stage$bottom_panel <- 
    clust_clinic$plots[c('cm10', 
                         'cm25ev', 
                         'cm25ni')] %>% 
    map2(., globals$cohort_lexicon$label[c('cm10', 
                                           'cm25ev', 
                                           'cm25ni')], 
         ~.x$MSKCC_risk + 
          theme(legend.position = 'none', 
                axis.title.y = element_blank()) + 
          labs(title = .y, 
               subtitle = .x$MSKCC_risk$labels$subtitle %>% 
                 stri_replace(regex = '^V.*,\\s{1}', 
                              replacement = ''))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv') %>% 
    plot_grid(get_legend(clust_clinic$plots$cm10$MSKCC_risk + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.8, 0.2))
  
  ## the entire figure
  
  suppl_figures$stage <- 
    plot_grid(suppl_figures$stage$upper_panel, 
              suppl_figures$stage$bottom_panel, 
              nrow = 2, 
              rel_heights = c(2, 1.25), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s11_cluster_stage_distribution', 
              ref_name = 'stage', 
              caption = paste('Distribution of tumor stages and MSKCC risk', 
                              'strata in the lymphokine clusters.'), 
              w = 180, 
              h = 210)
  
# Figure S12: checkmate response rate -------
  
  insert_msg('Figure S12: checkmate response rate')
  
  suppl_figures$response <- 
    clust_clinic$plots[c('cm10', 'cm25ev', 'cm25ni')] %>% 
    map2(., globals$cohort_lexicon$label[c('cm10', 
                                           'cm25ev', 
                                           'cm25ni')], 
         ~.x$biresponse + 
           labs(title = .y, 
                subtitle = .x$biresponse$labels$subtitle %>% 
                  stri_replace(regex = '^V.*,\\s{1}', 
                               replacement = '')) + 
           theme(legend.position = 'none', 
                 axis.title.x = element_blank())) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv') %>% 
    plot_grid(get_legend(clust_clinic$plots$cm10$biresponse + 
                           theme(legend.position = 'bottom')), 
              nrow = 2, 
              rel_heights = c(0.85, 0.15)) %>%
    as_figure(label = 'figure_s12_cluster_therapy_response', 
              ref_name = 'response', 
              caption = paste('Therapy response rates in the lymptokine', 
                              'clusters in the CheckMate cohorts.'), 
              w = 180, 
              h = 85)
  
# Figure S13: relapse-free survival in the lymphokine clusters ------
  
  insert_msg('Figure S13: relapse-free survival in the clusters')
  
  suppl_figures$clust_rfs <- clust_surv$km_plots$rfs %>% 
    map(~.x + 
          theme(legend.position = 'bottom') + 
          labs(color = NULL)) %>% 
    plot_grid(plotlist = ., 
              ncol = 2, 
              align = 'hv') %>% 
    as_figure(label = 'figure_s13_cluster_relapse_free_survival', 
              ref_name = 'clust_rfs', 
              caption = paste('Relapse-free survival in the', 
                              'lymphokine clusters.'), 
              w = 180, 
              h = 150)
  
# Figures S: univariable OS for the genes of interest -------
  
  ## not show in the paper but available in the repository
  ## for interested readers
  
  #insert_msg('Figures S14 - S15: univariable OS, genes of interest')
  
 # suppl_figures[c('os1', 'os2')] <- 
  #  list(uni_cox$forest[c('CXCL9', 'CXCL10', 'CXCL11', 'CXCR3')], 
   #      uni_cox$forest[c('CXCL13', 'CXCR5', 'XCL1', 'XCR1')]) %>% 
   # map(~map(.x, 
   #          ~.x + 
   #            theme(plot.subtitle = element_blank(), 
   #                  legend.position = 'none',
   #                  axis.text = element_text(size = 7)))) %>% 
  #  map(~plot_grid(plotlist = .x, 
   #                ncol = 2, 
    #               align = 'hv') %>% 
     #     plot_grid(get_legend(uni_cox$forest[[1]] + 
      #                           scale_color_manual(values = c(upregulated = 'firebrick', 
       #                                                        downregulated = 'steelblue', 
        #                                                       ns = 'gray60'), 
         #                                           labels = c(upregulated = 'unfavorable', 
          #                                                     downregulated = 'favorable', 
           #                                                    ns = 'ns'), 
            #                                        name = 'Effect on survival') + 
             #                    guides(fill = 'none') + 
              #                   theme(legend.position = 'bottom')), 
               #     nrow = 2, 
                #    rel_heights = c(0.93, 0.07)))
  
 # suppl_figures[c('os1', 'os2')] <- 
  #  suppl_figures[c('os1', 'os2')] %>% 
   # list(x = ., 
    #     label = c('figure_s14_overall_survival', 
     #              'figure_s15_overall_survival'), 
      #   ref_name = c('os1', 'os2'), 
       #  caption = c(paste('Modeling of overall survival as a', 
        #                   'function of cancer expression of CXCL9/10/11', 
         #                  'and CXCR3.'), 
          #           paste('Modeling of overall survival as a', 
           #                'function of cancer expression of CXCL13,', 
            #               'CXCR5, XCL1 and XCR1.'))) %>% 
  # pmap(as_figure, 
   #      w = 180, 
    #     h = 150)
  
# Figures S: univariable RFS for the genes of interest -------
  
  ## not show in the paper but available in the repository
  ## for interested readers
  
 # insert_msg('Figures S16 - S17: univariable RFS, genes of interest')
  
  #suppl_figures[c('rfs1', 'rfs2')] <- 
   # list(uni_rfs$forest[c('CXCL9', 'CXCL10', 'CXCL11', 'CXCR3')], 
    #     uni_rfs$forest[c('CXCL13', 'CXCR5', 'XCL1', 'XCR1')]) %>% 
  #  map(~map(.x, 
   #          ~.x + 
    #           theme(plot.subtitle = element_blank(), 
     #                legend.position = 'none',
      #               axis.text = element_text(size = 7)))) %>% 
  #  map(~plot_grid(plotlist = .x, 
   #                ncol = 2, 
    #               align = 'hv') %>% 
     #     plot_grid(get_legend(uni_cox$forest[[1]] + 
      #                           scale_color_manual(values = c(upregulated = 'firebrick', 
       #                                                        downregulated = 'steelblue', 
        #                                                       ns = 'gray60'), 
         #                                           labels = c(upregulated = 'unfavorable', 
          #                                                     downregulated = 'favorable', 
           #                                                    ns = 'ns'), 
            #                                        name = 'Effect on survival') + 
             #                    guides(fill = 'none') + 
              #                   theme(legend.position = 'bottom')), 
               #     nrow = 2, 
                #    rel_heights = c(0.93, 0.07)))
  
 # suppl_figures[c('rfs1', 'rfs2')] <- 
  #  suppl_figures[c('rfs1', 'rfs2')] %>% 
   # list(x = ., 
    #     label = c('figure_s16_relapse_free_survival', 
     #              'figure_s17_relapse_free_survival'), 
      #   ref_name = c('rfs1', 'rfs2'), 
       #  caption = c(paste('Modeling of relapse-free survival as a', 
        #                   'function of cancer expression of CXCL9/10/11', 
         #                  'and CXCR3.'), 
          #           paste('Modeling of relapse-free survival as a', 
           #                'function of cancer expression of CXCL13,', 
            #               'CXCR5, XCL1 and XCR1.'))) %>% 
  #  pmap(as_figure, 
   #      w = 180, 
    #     h = 130)
  
# Figure S14: multi-variable ridge Cox regression --------
  
  insert_msg('Figure S14: Multi-variable ridge Cox regression')
  
  suppl_figures$multi_os <- 
    plot_grid(clin_os$coef_plots$genes + 
                theme(legend.position = 'none', 
                      plot.subtitle = element_blank(), 
                      axis.text = element_text(size = 7)), 
              plot_grid(clin_os$c_index_plot + 
                          theme(legend.position = 'bottom', 
                                legend.text = element_text(size = 7), 
                                axis.text = element_text(size = 7),  
                                plot.subtitle = element_blank()), 
                        nrow = 2, 
                        rel_heights = c(0.7, 0.3)) + 
                labs(title = 'Predction quality'), 
              ncol = 2, 
              align = 'h',
              axis = 'tblr', 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s14_multi_paramater_overall_survival', 
              ref_name = 'multi_os', 
              caption = paste('Multi-parameter modeling of overall survival', 
                              'as a function of gene expression, demographic', 
                              'and staging data.'), 
              w = 180, 
              h = 140)
  
# Figure S15: signaling in the clusters ------
  
  insert_msg('Figure S15: signaling in the clusters')
  
  suppl_figures$signaling <- 
    plot_grid(dge_spia$bubble_plots$activated) %>% 
    as_figure(label = 'figure_s15_cluster_signaling_pathways', 
              ref_name = 'signaling', 
              caption = paste('Predicted activation of', 
                              'signaling pathways in the chemox high', 
                              'versus chemox low lymphokine cluster.'), 
              w = 180, 
              h = 140)
  
# Figures S16 - S17: interferon gamma signaling ------
  
  insert_msg('Figures S16 - S17: interferon gamma')
  
  suppl_figures[c('ifn1', 'ifn2')] <- 
    list(dge_ifn$volcano_plots[c('tcga', 'emtab1980', 
                                 'gse73731', 'gse167093')], 
         dge_ifn$volcano_plots[c('reca', 'cm10', 
                                 'cm25ev', 'cm25ni')]) %>% 
    map(~map(.x, ~.x + theme(legend.position = 'none'))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = 'hv') %>% 
          plot_grid(get_legend(dge_ifn$volcano_plots[[1]] + 
                                 theme(legend.position = 'bottom')), 
                    nrow = 2, 
                    rel_heights = c(0.93, 0.07)))
  
  suppl_figures[c('ifn1', 'ifn2')] <- 
    suppl_figures[c('ifn1', 'ifn2')] %>% 
    list(x = ., 
         label = c('figure_s16_interferon_gamma_signaling', 
                   'figure_s17_interferon_gamma_signaling'), 
         ref_name = c('ifn1', 'ifn2'), 
         caption = c(paste('Differential expression of interferon', 
                           'gamma signaling-associated genes between the', 
                           'lymphokine clusters in the TCGA KIRC,', 
                           'E-MTAB 1980, GSE73731 and GSE167093 cohorts.'), 
                     paste('Differential expression of interferon', 
                           'gamma signaling-associated genes between the', 
                           'lymphokine clusters in the CheckMate cohorts.'))) %>% 
    pmap(as_figure, 
         w = 180, 
         h = 180)
  
# Figures S18 - S19: exhaustion genes in the clusters -------
  
  insert_msg('Exhaustion genes in the clusters')
  
  suppl_figures[c('exhaustion1', 'exhaustion2')] <- 
    list(dge_exhaust$volcano_plots[c('tcga', 'emtab1980', 
                                     'gse73731', 'gse167093')], 
         dge_exhaust$volcano_plots[c('reca', 'cm10', 
                                     'cm25ev', 'cm25ni')]) %>% 
    map(~map(.x, ~.x + theme(legend.position = 'none'))) %>% 
    map(~plot_grid(plotlist = .x, 
                   ncol = 2, 
                   align = 'hv') %>% 
          plot_grid(get_legend(dge_exhaust$volcano_plots[[1]] + 
                                 theme(legend.position = 'bottom')), 
                    nrow = 2, 
                    rel_heights = c(0.93, 0.07)))
  
  suppl_figures[c('exhaustion1', 'exhaustion2')] <- 
    suppl_figures[c('exhaustion1', 'exhaustion2')] %>% 
    list(x = ., 
         label = c('figure_s18_exhaustion', 
                   'figure_s19_exhaustion'), 
         ref_name = c('exhaustion1', 'exhaustion2'), 
         caption = c(paste('Differential expression of T cell', 
                           'exhaustion-associated genes between the', 
                           'lymphokine clusters in the TCGA KIRC,', 
                           'E-MTAB 1980, GSE73731 and GSE167093 cohorts.'), 
                     paste('Differential expression of T cell', 
                           'exhaustion-associated genes between the', 
                           'lymphokine clusters in the CheckMate cohorts.'))) %>% 
    pmap(as_figure, 
         w = 180, 
         h = 180)
  
# Figure S20: checkpoint genes in the Checkmate cohorts -------
  
  insert_msg('Figure S20: checkpoint genes in the clusters')
  
  suppl_figures$checkpoint <- dge_exhaust$checkpoint_plots[c('cm10', 
                                                             'cm25ev', 
                                                             'cm25ni')] %>% 
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
    as_figure(label = 'figure_s20_clusters_checkpoint', 
              ref_name = 'checkpoint', 
              caption = paste('Differential expression of genes coding', 
                              'for clinically relevant immune checkpoint', 
                              'proteins in the lymphokine clusters', 
                              'in the CheckMate cohorts'), 
              w = 180, 
              h = 140)
  
# Figure S21: fatty acid oxidation and TCA ---------
  
  insert_msg('Figure S21: fatty acid oxidation')
  
  ## the total number of reaction will be provided in the plot
  ## legend
  
  suppl_figures$fao <- meta_details$fao_sign_plots %>% 
    map(~.x + 
          labs(subtitle = .x$labels$subtitle %>% 
                 stri_replace(regex = '^total.*,\\s{1}activated', 
                              replacement = 'activated'), 
               x = 'FAO reaction', 
               y = expression('log'[2] * ' fold regulation')) + 
          theme(legend.position = 'none')) %>% 
    c(list(get_legend(meta_details$fao_sign_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s21_fatty_acid_oxidation', 
              ref_name = 'fao', 
              caption = paste('Predicted differences in activity of', 
                              'fatty acid oxidation reactions between the', 
                              'lymphokine clusters.'), 
              w = 180, 
              h = 160)
  
# Figure S22: citric acid cycle --------
  
  insert_msg('Figure S22: TCA reaction activity')
  
  suppl_figures$tca <- meta_details$tca_detail_plots %>% 
    map(~.x + 
          scale_y_discrete(labels = tca_labeller) + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'none', 
                axis.text = element_text(size = 7)) + 
          labs(expression('log'[2] * ' fold regulation'))) %>% 
    c(list(get_legend(meta_details$tca_detail_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr')
  
  suppl_figures$tca <- 
    plot_grid(ggdraw() + 
                draw_image('./schemes/tca.png') + 
                theme(plot.margin = ggplot2::margin(t = 5, l = 10, 
                                                    r = 10, b = 5, 
                                                    unit = 'mm')), 
              suppl_figures$tca, 
              nrow = 2, 
              rel_heights = c(1, 4), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s22_citric_acid_cycle', 
              ref_name = 'tca', 
              caption = paste('Predicted differences in activity of', 
                              'citric acid cycle reactions between the', 
                              'lymphokine clusters.'), 
              w = 180, 
              h = 230)
  
# Figure S23: oxidative phosphorylation -------
  
 insert_msg('Figure S23: oxidative phosphorylation')
  
  suppl_figures$oxphos <- meta_details$oxphos_sign_plots %>% 
    map(~.x + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'none')) %>% 
    c(list(get_legend(meta_details$oxphos_sign_plots[[1]]))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') %>% 
    as_figure(label = 'figure_s23_oxidative_phosphorylation', 
              ref_name = 'oxphos', 
              caption = paste('Predicted differences in activity of', 
                              'oxidative phosphorylation complexes between the', 
                              'lymphokine clusters.'), 
              w = 180, 
              h = 160)
  
# Figure S24: tryptophan decay --------
  
  insert_msg('Figure S24: tryptophan decay')

  suppl_figures$trycats <- meta_details$trycats_detail_plots %>% 
    map(~.x + 
          scale_y_discrete(labels = trycats_labeller) + 
          theme(plot.subtitle = element_blank(), 
                legend.position = 'none', 
                axis.text = element_text(size = 7), 
                plot.margin = ggplot2::margin(t = 2, r = 2, 
                                              l = 2, b = 0, 
                                              unit = 'mm'))) %>% 
    c(list(get_legend(meta_details$trycats_detail_plots[[1]] + 
                        theme(legend.text = element_text(size = 7))))) %>% 
    plot_grid(plotlist = ., 
              ncol = 3, 
              align = 'hv', 
              axis = 'tblr') %>% 
    plot_grid(ggdraw() + 
                draw_image('./schemes/trycats.png') + 
                theme(plot.margin = ggplot2::margin(t = 5, l = 10, 
                                                    r = 10, b = 5, 
                                                    unit = 'mm')), 
              ., 
              nrow = 2, 
              rel_heights = c(1, 4), 
              labels = LETTERS, 
              label_size = 10) %>% 
    as_figure(label = 'figure_s24_tryptophan_decay', 
              ref_name = 'trycats', 
              caption = paste('Predicted differences in activity of', 
                              'the tryptophan/kynurenine/quinolinate pathway', 
                              'between the lymphokine clusters.'), 
              w = 180, 
              h = 210)
  
# Figure S25: TRYCATS genes in the clusters, CheckMate cohorts ------
  
  insert_msg('Figure S25: TRYCATS genes in the clusters, CheckMate')
  
  suppl_figures$trycats_genes <- meta_genes$trycats_plots[c('cm10', 
                                                            'cm25ev', 
                                                            'cm25ni')] %>% 
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
    as_figure(label = 'figure_s25_trycats_genes', 
              ref_name = 'trycats_genes', 
              caption = paste('Differential expression of genes involved', 
                              'in the tryptophan/kynurenine/quinolinate', 
                              'metabolism in the lymphokine clusters in', 
                              'the CheckMate cohorts.'), 
              w = 180, 
              h = 140)
  
# saving the figures ------
  
  insert_msg('Saving the figures')
  
  suppl_figures %>% 
    walk(pickle, 
         path = './paper/supplementary figures', 
         format = 'pdf', 
         device = cairo_pdf)
  
# END ------
  
  insert_tail()