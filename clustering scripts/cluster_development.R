# Finding of the optimal clustering algorithm in the TCGA training cohort
#
# The criteria: optimal clustering variance and accuracy 
# in 10-fold cross-validation

  insert_head()
  
# container -------
  
  clust_dev <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  clust_dev$genes <- globals$cxcl_genes
  
  ## analysis tables: mean-centered gene expression
  
  clust_dev$analysis_tbl <- list()
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('clust_dev$analysis_tbl$', i, 
                           ' <- ', i, '$expression')))

    
  }
  
  clust_dev$analysis_tbl[c('tcga', 'gse167093')] <- 
    clust_dev$analysis_tbl[c('tcga', 'gse167093')] %>% 
    map(filter, tissue_type == 'Tumor')
  
  clust_dev$analysis_tbl <- clust_dev$analysis_tbl %>% 
    map(select, 
        patient_id, 
        all_of(clust_dev$genes)) %>% 
    map(column_to_rownames, 'patient_id') %>% 
    map(center_data, 'mean')
  
  ## distances to be tested
  
  clust_dev$distances <- c('euclidean', 'manhattan', 'cosine')
  
# algorithms -------
  
  insert_msg('Creating clustering objects')
  
  ## hierarchical clustering, the number of clusters based 
  ## on the WSS curve and the dendrogram
  
  clust_dev$algos[c(paste0('hcl_', clust_dev$distances))] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 3, 2)) %>% 
    pmap(hcluster, 
         data = clust_dev$analysis_tbl$tcga, 
         hc_method = 'ward.D2')
  
  ## Kmeans and PAM clustering
  
  clust_dev$algos[c(paste0('kmeans_', clust_dev$distances))] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$tcga, 
         clust_fun = 'kmeans')
  
  clust_dev$algos[c(paste0('pam_', clust_dev$distances))] <- 
    list(distance_method = clust_dev$distances, 
         k = c(3, 3, 2)) %>% 
    pmap(kcluster, 
         data = clust_dev$analysis_tbl$tcga, 
         clust_fun = 'pam')
  
# clustering variance -----
  
  insert_msg('Clustering variance')
  
  clust_dev$variance <- clust_dev$algos %>% 
    map(var) %>% 
    map_dbl(~.x$frac_var) %>% 
    compress(names_to = 'method', 
             values_to = 'variance')
  
# cross-validation -------
  
  insert_msg('Cross validation')
  
  ## 10-fold cross-validation with re-assignment
  ## with inverse distance-weighted 7-NN classifier

  clust_dev$cv <- clust_dev$algos %>% 
    map(cv, 
        nfolds = 10, 
        kNN = 7, 
        simple_vote = FALSE, 
        resolve_ties = TRUE, 
        .parallel = TRUE)
  
  clust_dev$cv <- clust_dev$cv %>% 
    map(~.x$summary) %>% 
    compress(names_to = 'method') %>% 
    mutate(accuracy = 1 - mean_error)
  
# Testing results -------
  
  insert_msg('Testing results')
  
  clust_dev$test_results <- 
    left_join(clust_dev$variance, 
              clust_dev$cv, 
              by = 'method') %>% 
    mutate(algorithm = stri_split_fixed(method, 
                                        pattern = '_', 
                                        simplify = TRUE)[, 1], 
           algorithm = toupper(algorithm), 
           distance = stri_split_fixed(method, 
                                       pattern = '_', 
                                       simplify = TRUE)[, 2], 
           method_lab = paste(algorithm, distance, sep = ', '))
  
# Plotting the testing results -------
  
  insert_msg('Plotting the testing results')
  
  clust_dev$result_plot <- clust_dev$test_results %>% 
    pivot_longer(cols = c('variance', 'accuracy'), 
                 names_to = 'stat', 
                 values_to = 'value') %>% 
    ggplot(aes(x = value, 
               y = reorder(method_lab, value), 
               fill = stat)) + 
    geom_bar(stat = 'identity', 
             color = 'black', 
             position = position_dodge(0.9)) + 
    geom_text(aes(label = signif(value, 2)), 
              size = 2.6, 
              color = 'white', 
              hjust = 1.4, 
              position = position_dodge(0.9)) + 
    scale_fill_manual(values = c(accuracy = 'darkolivegreen4', 
                                 variance = 'steelblue3'), 
                      labels = c(accuracy = 'CV accuracy', 
                                 variance = 'clustering variance'), 
                      name = '') + 
    globals$common_theme + 
    theme(axis.title.y = element_blank()) + 
    labs(title = 'Clustering by lymphokine genes, TCGA KIRC', 
         subtitle = 'Explained clusteing variance and 10-fold cross-validation', 
         x = 'statistic value')
  
# saving the trained clusters ------
  
  insert_msg('Saving the trained clusters')
  
  save(clust_dev, file = './cache/gene_clusters.RData')
  
# END -----
  
  rm(i)
  
  insert_tail()