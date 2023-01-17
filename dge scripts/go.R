# GO enrichment for the differentially regulated genes

  insert_head()
  
# container -----
  
  dge_go <- list()
  
# parallel backend -----
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# enrichment ------
  
  insert_msg('Enrichment')
  
  dge_go$test_results <- clust_dge$entrez_vectors %>% 
    map(unname) %>% 
    future_map(~goana(de = .x), 
               .options = furrr_options(seed = TRUE)) %>% 
    map(rownames_to_column, 'go_id') %>% 
    map(filter, Ont == 'BP') %>% 
    map(mutate, p_adjusted = p.adjust(P.DE, 'BH')) %>% 
    map(filter, P.DE < 0.05) %>% 
    map(as_tibble)
  
# common regulated GOs -------
  
  insert_msg('Common significantly enriched GOs')
  
  dge_go$common_go <- dge_go$test_results %>% 
    map(filter, p_adjusted < 0.05) %>% 
    map(~.x$Term) %>% 
    reduce(intersect)
  
# Top enriched GOs -------
  
  insert_msg('Top enriched GOs')
  
  ## within the set of commonly enriched GOs
  
  dge_go$enrich_plots <- 
    list(data = dge_go$test_results %>% 
           map(filter, Term %in% dge_go$common_go), 
         plot_title = globals$cohort_lexicon$label) %>% 
    pmap(plot_signifcant, 
         p_variable = 'p_adjusted', 
         label_variable = 'Term', 
         signif_level = 0.05, 
         top_significant = 20, 
         x_lab = expression('-log'[10] * ' pFDR'), 
         globals$common_theme)
  
# END ------
  
  plan('sequential')
  
  insert_tail()