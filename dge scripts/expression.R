# Identification of genes differentially expressed between the lymphokine 
# clusters
#
# Differentially regulated genes are identified
# by two-tailed FDR-corrected T test and at least 1.5-fold regulation between 
# the clusters

  insert_head()
  
# container -------
  
  clust_dge <- list()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## expression tables and analysis genes
  
  clust_dge$annotation <- list()
  clust_dge$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('clust_dge$annotation$', i, 
                           ' <- ', i, '$annotation')))
    
    eval(parse_expr(paste0('clust_dge$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))
    
    
  }
  
  clust_dge$annotation <- clust_dge$annotation %>% 
    map(~filter(.x, 
                complete.cases(.x), 
                gene_symbol != '', 
                entrez_id != '')) %>% 
    map2(., clust_dge$analysis_tbl, 
         ~filter(.x, gene_symbol %in% names(.y))) %>% 
    map(filter, !duplicated(gene_symbol))
  
  clust_dge$variables <- clust_dge$annotation %>% 
    map(~.x$gene_symbol)
  
  clust_dge$analysis_tbl[c('tcga', 'gse167093')] <- 
    clust_dge$analysis_tbl[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  clust_dge$analysis_tbl <- 
    map2(clust_dge$analysis_tbl, 
         clust_dge$variables, 
         ~.x[c('patient_id', .y)])
  
  ## appending the expression tables with the clustering information
  
  clust_dge$analysis_tbl <- gene_clust$clust_obj %>% 
    map(~.x$clust_assignment) %>% 
    map(set_names, c('patient_id', 'clust_id')) %>% 
    future_map2(., clust_dge$analysis_tbl, 
                left_join, by = 'patient_id')
  
# serial testing ---------
  
  insert_msg('Serial testing')
  
  clust_dge$test_results <- 
    list(data = clust_dge$analysis_tbl, 
         variables = clust_dge$variables) %>% 
    pmap(test_two_groups, 
         split_fct = 'clust_id', 
         adj_method = 'BH', 
         .parallel = TRUE)
  
  ## annotation of the results with Entrez IDs
  
  clust_dge$test_results <- 
    map2(clust_dge$test_results, 
         clust_dge$annotation, 
         ~mutate(.x, 
                 gene_symbol = response, 
                 entrez_id = exchange(gene_symbol, 
                                      dict = .y, 
                                      key = 'gene_symbol', 
                                      value = 'entrez_id'), 
                 regulation = ifelse(significant == 'no', 
                                     'ns', 
                                     ifelse(estimate > log2(1.25), 
                                            'upregulated', 
                                            ifelse(estimate < -log2(1.25), 
                                                   'downregulated', 'ns'))), 
                 regulation = factor(regulation, 
                                     c('upregulated', 'ns', 'downregulated'))))

# Identification of differentially regulated genes -------
  
  insert_msg('Identification of differentially regulated genes')
  
  clust_dge$significant <- clust_dge$test_results %>% 
    map(filter, regulation != 'ns')
  
  ## common up- and downregulated genes 
  ## (all cohorts except of the smallest CM10)
  
  clust_dge$common_up <- 
    clust_dge$significant[names(clust_dge$significant) != 'cm10'] %>% 
    map(filter, regulation == 'upregulated') %>% 
    map(~.x$gene_symbol) %>% 
    reduce(intersect)
  
  clust_dge$common_down <- 
    clust_dge$significant[names(clust_dge$significant) != 'cm10'] %>% 
    map(filter, regulation == 'downregulated') %>% 
    map(~.x$gene_symbol) %>% 
    reduce(intersect)
  
# Entrez ID lists and fold-regulation vectors ------
  
  insert_msg('Enrez ID lists and fold-regulation vectors')
  
  ## lists of significantly regulated Entrez IDs
  
  clust_dge$entrez_vectors <- clust_dge$significant %>% 
    map(~.x$entrez_id)
  
  clust_dge$regulation_vectors <- clust_dge$significant %>% 
    map(~set_names(.x$estimate, 
                   .x$entrez_id))
  
# saving the testing results ------
  
  insert_msg('Saving the testing results')
  
  save(clust_dge, file = './cache/dge.RData')
  
# END -------
  
  rm(i)
  
  plan('sequential')
  
  insert_tail()