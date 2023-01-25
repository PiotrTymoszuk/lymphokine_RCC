# Multi-dimensional scaling to identify the nearest partnets of CXCL9, CXCL10, 
# CXCL11 among xCell immune infiltration estimates

  insert_head()

# container ------

  xcell_mds <- list()

# analysis globals --------

  insert_msg('Analysis globals')

  ## analysis tables: median centering

  xcell_mds$analysis_tbl <- list()
  xcell_mds$expression_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('xcell_mds$analysis_tbl$', i, 
                           ' <- ', i, '$xcell')))
    
    eval(parse_expr(paste0('xcell_mds$expression_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  xcell_mds$analysis_tbl <- xcell_mds$analysis_tbl %>% 
    map(~.x[c('patient_id', globals$xcell_lexicon$variable)])
  
  xcell_mds$expression_tbl <- xcell_mds$expression_tbl %>% 
    map(select, any_of(c('patient_id', 'tissue_type', globals$cxcl_genes)))
  
  xcell_mds$expression_tbl[c('tcga', 'gse167093')] <- 
    xcell_mds$expression_tbl[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  xcell_mds$expression_tbl <- xcell_mds$expression_tbl %>% 
    map(select, patient_id, all_of(globals$cxcl_genes))
  
  xcell_mds$analysis_tbl <- 
    map2(xcell_mds$analysis_tbl, 
         xcell_mds$expression_tbl, 
         inner_join, by = 'patient_id') %>% 
    map(column_to_rownames, 'patient_id') %>% 
    #map(center_data, 'median') %>% 
    map(as.data.frame)
  
# Multi-dimensional scaling: median-centering, distance calculation ------
  
  insert_msg('Multi-dimensional scaling')
  
  ## two-dimensional MDS with Spearman distance between the observations
  ## calculating by hand since the  scaled Spearman distance is not implemented
  ## in the clustTool package
  
  ## distance objects
  
  xcell_mds$dist_mtx <- xcell_mds$analysis_tbl %>% 
    map(cor, method = 'spearman') %>% 
    map(function(x) max(x) - x) %>% 
    map(as.dist, diag = TRUE)

  xcell_mds$mds_obj <- xcell_mds$dist_mtx %>% 
    #map(as.matrix) %>% 
    #map(umap::umap, input = 'dist', n_components = 3) %>% 
    #map(~.x$layout) %>% 
    map(cmdscale, k = 3) %>% 
    map(as.data.frame) %>% 
    map(~set_names(.x, 
                   paste0('comp_', 1: ncol(.x)))) %>% 
    map(rownames_to_column, 'observation') %>% 
    as_tibble
  
  xcell_mds$mds_obj <- 
    list(x = map(xcell_mds$analysis_tbl, ~quo(.x)), 
         y = xcell_mds$mds_obj) %>% 
    pmap(function(x, y) list(red_obj = NULL, 
                             red_fun = 'mds', 
                             component_tbl = y, 
                             loadings = NULL, 
                             data = x)) %>% 
    map(red_analysis)

# Identification of the nearest neighbors of the genes of interest -------
  
  insert_msg('Nearest neighbors')

  ## nearest neighbors: selecting the cell populations
  ## associated with the chemokine gene with rho > 0.4
  
  xcell_mds$knn_labels <- xcell_mds$dist_mtx %>% 
    map(as.matrix) %>% 
    map(~.x[!rownames(.x) %in% globals$cxcl_genes, globals$cxcl_genes]) %>% 
    map(as.data.frame) %>% 
    map(rownames_to_column, 'cell_type') %>% 
    map(function(cohort) set_names(globals$cxcl_genes, 
                                   globals$cxcl_genes) %>% 
          map(~filter(cohort[, c('cell_type', .x)], .data[[.x]] < 0.6)) %>% 
          map(~.x$cell_type)) %>% 
    map(compress, 
        names_to = 'gene', 
        values_to = 'knn')

# Plotting the MDS layouts with the nearest neighbors labeled ------
  
  insert_msg('Plotting')

  xcell_mds$mds_plots <- list(mds_object = xcell_mds$mds_obj, 
                              knn_data = xcell_mds$knn_labels, 
                              plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_mds, 
         nearest_only = FALSE)
  
# Force-directed graph of scaled Spearman's correlation coefficients ----
  
  insert_msg('Force-directed graphs')
  
  set.seed(12345)
  
  xcell_mds$fd_plots <- 
    list(data = xcell_mds$analysis_tbl, 
         cell_subset = xcell_mds$knn_labels %>% 
           map(~.x$knn) %>% 
           map(reduce, union), 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(forced_netgraph, 
         weighting_order = 6, 
         gene_color = 'firebrick4', 
         cell_color = 'steelblue4')

# END ------
  
  insert_tail()