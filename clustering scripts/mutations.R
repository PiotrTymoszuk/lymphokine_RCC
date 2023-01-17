# Checking for differences in somatic mutation burden and frequency between
# the lymphokine clusters. The mutation data are only available for
# the TCGA cohort.

  insert_head()
  
# container ------
  
  clust_mutations <- list()
  
# analysis globals ---------
  
  insert_msg('Analysis globals')
  
  ## total mutation burden
  
  clust_mutations$analysis_tbl$tmb <- 
    inner_join(somut$tmb, 
               gene_clust$clust_obj$tcga$clust_assignment %>% 
                 set_names(c('patient_id', 'clust_id')), 
               by = 'patient_id')
  
  ## frequencies of mutated alleles

  clust_mutations$analysis_tbl$freq <- 
    inner_join(somut$mut_summary, 
               gene_clust$clust_obj$tcga$clust_assignment %>% 
                 set_names(c('patient_id', 'clust_id')), 
               by = 'patient_id') %>% 
    map_dfc(function(x) if(is.numeric(x)) factor(car::recode(x, 
                                                      "'0' = 'WT'; '1' = 'mutated'"), 
                                                 c('WT', 'mutated')) else x)
  
  ## frequent mutations to be analyzed 
  
  clust_mutations$variables <- somut$top_mutations$variable
  
# TMB ------
  
  insert_msg('TMB')
  
  ## descriptive stats
  
  clust_mutations$tmb_stats <- 
    clust_mutations$analysis_tbl$tmb %>% 
    explore(split_factor = 'clust_id', 
            variables = 'mut_count', 
            what = 'table', 
            pub_styled = TRUE)
  
  ## Wilcoxon test
  
  clust_mutations$tmb_test <- 
    clust_mutations$analysis_tbl$tmb %>% 
    compare_variables(split_factor = 'clust_id', 
                      variables = 'mut_count', 
                      what = 'eff_size',
                      types = 'wilcoxon_r', 
                      ci = FALSE, 
                      exact = FALSE, 
                      pub_styled = TRUE) %>% 
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))
  
  ## plotting 
  
  clust_mutations$tmb_plot <- 
    plot_variable(clust_mutations$analysis_tbl$tmb, 
                  variable = 'mut_count', 
                  split_factor = 'clust_id', 
                  type = 'violin', 
                  point_hjitter = 0, 
                  cust_theme = globals$common_theme, 
                  plot_title = 'TMB', 
                  plot_subtitle = clust_mutations$tmb_stats$plot_cap, 
                  x_lab = 'Lymphokine cluster', 
                  y_lab = '# mutated genes') + 
    scale_fill_manual(values = globals$cluster_colors) + 
    guides(fill = 'none') + 
    scale_y_continuous(trans = 'log2')
  
# Differences in mutation frequency --------
  
  insert_msg('Differences in mutation frequency')
  
  ## descriptive stats
  
  clust_mutations$freq_stats <- clust_mutations$analysis_tbl$freq %>% 
    explore(split_factor = 'clust_id', 
            variables = clust_mutations$variables, 
            what = 'table', 
            pub_styled = TRUE)
  
  ## testing results 
  
  clust_mutations$freq_test <- clust_mutations$analysis_tbl$freq %>% 
    compare_variables(split_factor = 'clust_id', 
                      variables = clust_mutations$variables, 
                      what = 'eff_size',
                      types = 'cramer_v',
                      ci = FALSE, 
                      exact = FALSE, 
                      pub_styled = TRUE, 
                      adj_method = 'BH')
  
  ## near-significant genes
  
  clust_mutations$freq_significant <- clust_mutations$freq_test %>% 
    filter(p_value < 0.05)
  
  ## plots
  
  clust_mutations$freq_significant_plots <- 
    list(variable = clust_mutations$freq_significant$variable, 
         plot_title = exchange(clust_mutations$freq_significant$variable, 
                               dict = somut$top_mutations, 
                               key = 'variable', 
                               value = 'gene'), 
         plot_subtitle = clust_mutations$freq_significant$significance) %>% 
    pmap(plot_variable, 
         clust_mutations$analysis_tbl$freq, 
         split_factor = 'clust_id', 
         scale = 'percent', 
         type = 'stack', 
         cust_theme = globals$common_theme, 
         x_lab = 'Lymphokine cluster', 
         y_lab = '% of cluster') %>% 
    map(~.x + 
          scale_fill_brewer() + 
          theme(plot.title = element_text(face = 'italic'))) %>% 
    set_names(clust_mutations$freq_significant$variable)
    
# END -------
  
  insert_tail()