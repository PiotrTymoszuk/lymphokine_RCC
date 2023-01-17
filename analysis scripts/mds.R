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
    map(center_data, 'median') %>% 
    map(t) %>% 
    map(as.data.frame)
  
# Multi-dimensional scaling: median-centering, distance calculation ------
  
  insert_msg('Multi-dimensional scaling')
  
  ## two-dimensional MDS with cosine distance between the immune/gene features
  
  xcell_mds$mds_obj <- xcell_mds$analysis_tbl %>% 
    map(reduce_data, 
        distance_method = 'cosine', 
        kdim = 2, 
        red_fun = 'mds')
  
# Identification of the nearest neighbors of the genes of interest -------
  
  insert_msg('Nearest neighbors')
  
  ## distance matrix and search for the nearest neighbors
  
  xcell_mds$dist_mtx <- xcell_mds$analysis_tbl %>% 
  map(calculate_dist, 
        method = 'cosine') %>% 
    map(as.dist)
  
  ## nearest neighbor objects
  
  xcell_mds$knn <- xcell_mds$dist_mtx %>% 
    map(kNN, k = 10)
  
  ## extracting names (labels) of the 5 nearest immune 
  ## neighbors of the genes of interest
  
  xcell_mds$knn_labels <- xcell_mds$knn %>% 
    map(~.x$id) %>% 
    map(rownames)
  
  xcell_mds$knn_labels <- xcell_mds$knn %>% 
    map(~.x$id[globals$cxcl_genes, ]) %>% 
    map2(., xcell_mds$knn_labels, 
         ~matrix(.y[.x], 
                 nrow = length(globals$cxcl_genes), 
                 ncol = 10, 
                 dimnames = list(globals$cxcl_genes, 1:10)))
  
  xcell_mds$knn_labels <- xcell_mds$knn_labels %>% 
    map(function(cohort) globals$cxcl_genes %>% 
          map(~cohort[.x, ]) %>% 
          map(~.x[!.x %in% globals$cxcl_genes]) %>% 
          map(~.x[1:5]) %>% 
          map(unname) %>% 
          set_names(globals$cxcl_genes)) %>% 
    map(compress, 
        names_to = 'gene', 
        values_to = 'knn')
  
# Plotting the MDS layouts with the nearest neighbors labeled ------
  
  insert_msg('Plotting')

  xcell_mds$mds_plots <- list(mds_object = xcell_mds$mds_obj, 
                              knn_data = xcell_mds$knn_labels, 
                              plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_mds)

# END ------
  
  insert_tail()