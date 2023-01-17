# UMAP xCell infiltration estimates, cosine distance
# Overlay with expression of the genes of interest

  insert_head()
  
# container ------
  
  xcell_clust <- list()
  
# analysis globals --------

  insert_msg('Analysis globals')
  
  ## analysis tables: median centering
  
  xcell_clust$analysis_tbl <- list()
  xcell_clust$expression_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('xcell_clust$analysis_tbl$', i, 
                           ' <- ', i, '$xcell')))
    
    eval(parse_expr(paste0('xcell_clust$expression_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  xcell_clust$analysis_tbl <- xcell_clust$analysis_tbl %>% 
    map(~.x[c('patient_id', globals$xcell_lexicon$variable)])
  
  xcell_clust$expression_tbl <- xcell_clust$expression_tbl %>% 
    map(select, any_of(c('patient_id', 'tissue_type', globals$cxcl_genes)))
  
  xcell_clust$expression_tbl[c('tcga', 'gse167093')] <- 
    xcell_clust$expression_tbl[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
# construction of UMAP objects ------
  
  insert_msg('Construction of UMAP objects')
  
  plan('multisession')
  
  xcell_clust$umap_obj <- xcell_clust$analysis_tbl %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(center_data, 'median') %>% 
    future_map(~reduce_data(data = .x, 
                     distance_method = 'cosine', 
                     kdim = 2, 
                     red_fun = 'umap', 
                     nepochs = 200, 
                     min_dist = 0.1, 
                     random_state = 1234), 
               .options = furrr_options(seed = TRUE))
  
  plan('sequential')
  
# Extraction of the layouts and appending them with expression values ------
  
  insert_msg('UMAP layouts')
  
  xcell_clust$layouts <- xcell_clust$umap_obj %>% 
    map(extract, 'scores') %>% 
    map(mutate, patient_id = observation)

  xcell_clust$layouts <- 
    map2(xcell_clust$layouts, 
         xcell_clust$analysis_tbl, 
         left_join, by = 'patient_id')
  
  xcell_clust$layouts <- 
    map2(xcell_clust$layouts, 
         xcell_clust$expression_tbl, 
         left_join, by = 'patient_id')
  
# Plots for CXCL9, 10, 11 and CXCR3 -------
  
  insert_msg('UMAP layout plots')
  
  ## plot titles
  
  xcell_clust$gene_titles <- globals$cohort_lexicon$label %>%
    map(function(cohort) globals$cxcl_genes %>% 
          map(~paste0('<b><em>', .x, '</em>, ', cohort, '</b>')))
  
  ## UMAP plots
  
  xcell_clust$gene_layout_plots <- list(cohort = xcell_clust$layouts, 
                                        title = xcell_clust$gene_titles) %>% 
    pmap(function(cohort, title) list(overlay = globals$cxcl_genes, 
                                      plot_title = title) %>% 
           pmap(plot_umap, 
                data = cohort, 
                center = 'median', 
                point_circle = 'none', 
                point_size = 1) %>% 
           map(~.x + theme(plot.title = element_markdown())) %>% 
           set_names(globals$cxcl_genes))
  
# Plots of the immune infiltrates --------
  
  insert_msg('Plots for immune infiltrates')
  
  ## plot titles
  
  xcell_clust$immune_titles <- globals$cohort_lexicon$label %>%
    map(function(cohort) globals$xcell_lexicon$label %>% 
          map(~paste(.x, cohort, sep = ', ')))
  
  ## UMAP plots
  
  xcell_clust$immune_layout_plots <- list(cohort = xcell_clust$layouts, 
                                          title = xcell_clust$immune_titles) %>% 
    pmap(function(cohort, title) list(overlay = globals$xcell_lexicon$variable, 
                                      plot_title = title) %>% 
           pmap(plot_umap, 
                data = cohort, 
                center = 'median', 
                point_circle = 'none', 
                point_size = 1) %>% 
           set_names(globals$xcell_lexicon$variable))

# END ------
  
  rm(i)
  
  insert_tail()