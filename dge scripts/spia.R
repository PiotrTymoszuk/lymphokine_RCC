# Signaling pathway analysis with SPIA

  insert_head()
  
# container ------
  
  dge_spia <- list()
  
# parallel backend -------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# analysis globals --------
  
  insert_msg('Analysis globals')
  
  ## there are some duplicated entries in the GSE167093
  ## regulated gene vector
  
  dge_spia$de_vectors <- clust_dge$regulation_vectors %>% 
    map(~.x[!duplicated(names(.x))])
  
  dge_spia$all_vectors <- clust_dge$test_results %>% 
    map(~.x$entrez_id) %>% 
    map(unique)

# Enrichment analysis -----
  
  insert_msg('Enrichment analysis')
  
  dge_spia$test_results <- 
    list(de = dge_spia$de_vectors, 
         all = dge_spia$all_vectors) %>% 
    future_pmap(spia, 
                verbose = FALSE, 
                .options = furrr_options(seed = TRUE))
  
  dge_spia$test_results <- dge_spia$test_results %>% 
    map(as_tibble)
  
# Common regulated pathways --------
  
  insert_msg('Common regulated pathways')
  
  ## for all cohorts except for the smallest CM10
  
  dge_spia$common_pathways$activated <- 
    dge_spia$test_results[names(dge_spia$test_results) != 'cm10'] %>% 
    map(filter, 
        pGFdr < 0.05, 
        tA > 0) %>% 
    map(~.x$Name) %>% 
    reduce(intersect)
  
  dge_spia$common_pathways$inhibited <- 
    dge_spia$test_results[names(dge_spia$test_results) != 'cm10'] %>% 
    map(filter, 
        pGFdr < 0.05, 
        tA < 0) %>% 
    map(~.x$Name) %>% 
    reduce(intersect)
  
# Regulation plots with the regulation estimates --------
  
  insert_msg('Regulation plots for the common regulated pathways')
  
  ## no common activated pathways were identified
  
  dge_spia$regulation_plots$activated <- 
    list(data = dge_spia$test_results %>% 
           map(filter, Name %in% dge_spia$common_pathways$activated), 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_top, 
         regulation_variable = 'tA', 
         label_variable = 'Name', 
         p_variable = 'pGFdr', 
         signif_level = 0.05, 
         regulation_level = 0,
         top_regulated = 100, 
         x_lab = 'Regulation, chemox high vs low', 
         cust_theme = globals$common_theme)
  
# Bubble plot with the regulation estimates --------
  
  insert_msg('Regulation estimates')
  
  ## bubble plot for the activated pathways only
  ## no common inhibited pathways were identified
  
  dge_spia$bubble_plots$activated <- dge_spia$test_results %>% 
    map(filter, Name %in% dge_spia$common_pathways$activated) %>% 
    map(select, Name, tA, pGFdr) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(significant = ifelse(pGFdr < 0.05, 'p < 0.05', 'ns'), 
           cohort = factor(cohort, globals$cohort_lexicon$cohort), 
           tA = ifelse(significant != 'ns', tA, NA)) %>% 
    ggplot(aes(x = cohort, 
               y = reorder(Name, tA, na.rm = TRUE), 
               fill = tA, 
               size = abs(tA))) + 
    geom_point(shape = 21) + 
    geom_text(aes(label = signif(tA, 2)), 
               size = 2.4, 
               hjust = 0.5, 
               vjust = -1.5) + 
    scale_fill_gradient2(low = 'steelblue', 
                         mid = 'white', 
                         high = 'firebrick', 
                         midpoint = 0, 
                         name = 'Regulation\nchemox high vs low') + 
    scale_size_area(name = 'Regulation\nchemox high vs low') + 
    scale_x_discrete(labels = globals$cohort_lexicon$label) + 
    globals$common_theme + 
    theme(axis.title = element_blank(), 
          axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
    labs(title = 'Common activated signaling pathways', 
         subtitle = 'chemox high vs chemox low cancers')
  
# saving the results -------
  
  insert_msg('Saving the results')
  
  save(dge_spia, file = './cache/spia.RData')

# END ------
  
  plan('sequential')
  
  insert_tail()