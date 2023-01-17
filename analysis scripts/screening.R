# Screening for correlation of chemokine gene expression with the CD8+ T cell
# content predicted by QuanTIseq. Comparison of chemokine gene expression 
# between RCC and the normal kidney to identify genes specific for malignancy

  insert_head()
  
# container --------
  
  screen <- list()
  
# analysis globals ------
  
  insert_msg('Analysis globals')
  
  ## Chemokine genes
  ## names retrieved from HGNC 
  ## (https://www.genenames.org/data/genegroup/#!/group/483)
  ## restricting to those present in the TCGA dataset
  ## and the remaining datasets
  
  screen$variables <- read_csv('./data/chemokine_genes.csv', skip = 1)
  
  screen$variables <- 
    screen$variables[c('Approved symbol', 'Previous symbols')] %>% 
    map(stri_split, fixed = ',') %>% 
    map(unlist) %>% 
    map(~.x[!is.na(.x)]) %>% 
    reduce(union)
  
  for(i in globals$cohort_lexicon$cohort) {
    
    eval(parse_expr(paste0('screen$variables <-', 
                           'screen$variables[screen$variables %in% names(', 
                           i, 
                           '$expression)]')))
    
    
  }
  
  ## analysis tables with the chemokine gene expression values
  ## and quantiseq estimates
  
  screen$analysis_tbl$tissue <- tcga$expression[c('patient_id', 
                                                  'tissue_type', 
                                                  screen$variables)]
  
  screen$analysis_tbl$t_cells <- tcga$quantiseq[c('patient_id', 
                                                  'T cell CD8+')] %>% 
    left_join(tcga$expression[c('patient_id', 
                                'tissue_type', 
                                screen$variables)], 
              by = 'patient_id')
  
  screen$analysis_tbl$t_cells <- 
    screen$analysis_tbl$t_cells %>% 
    filter(tissue_type == 'Tumor')
  
  ## restricting the tissue expression table to the paired samples
  
  screen$tissue_pair_ids <- screen$analysis_tbl$tissue %>% 
    filter(tissue_type == 'Normal') %>% 
    .$patient_id
  
  screen$analysis_tbl$tissue <- screen$analysis_tbl$tissue %>% 
    filter(patient_id %in% screen$tissue_pair_ids) %>% 
    arrange(patient_id)
  
# Comparison of expression between the tumor and normal kidney -------
  
  insert_msg('Normal - tumor comparison')
  
  ## testing: two-tailed T test with the log2 transformed expression
  
  screen$norm_tumor_test <- 
    test_two_groups(screen$analysis_tbl$tissue, 
                    split_fct = 'tissue_type', 
                    variables = screen$variables, 
                    type = 't', 
                    adj_method = 'BH', 
                    paired = TRUE) %>% 
    mutate(estimate = -expect_x, 
           expect_y = 0)
  
  ## n numbers to be presented in the plot tag
  
  screen$norm_tumor_n <- 
    paste('\ntissue pairs: n =', length(screen$tissue_pair_ids))
  
  ## top genes
  
  screen$norm_tumor_top <- screen$norm_tumor_test %>% 
    group_by(sign(estimate)) %>% 
    top_n(n = 10, abs(estimate)) %>% 
    ungroup %>% 
    filter(significant == 'yes') %>% 
    mutate(fontface = ifelse(response %in% globals$cxcl_genes, 
                             'bold.italic', 
                             'italic'), 
           alpha = ifelse(response %in% globals$cxcl_genes, 
                          'high', 'low'), 
           regulation = ifelse(estimate > 0, 'upregulated', 'downregulated'))
  
  ## volcano plot
  
  screen$norm_tumor_volcano <- 
    plot_volcano(screen$norm_tumor_test, 
                 regulation_variable = 'estimate', 
                 p_variable = 'p_adjusted', 
                 signif_level = 0.05, 
                 regulation_level = 0, 
                 x_lab = expression('log'[2] * ' regulation, RCC vs normal kidney'), 
                 y_lab = expression('-log'[10] * 'pFDR'), 
                 cust_theme = globals$common_theme, 
                 plot_title = paste('Chemokine genes, RCC vs normal kidney,', 
                                    'TCGA KIRC'))
  
  ## appending with thr total gene numbers and labeling the top genes
  ## CXCL9/10/11 in bold/italic
  
  screen$norm_tumor_volcano <- screen$norm_tumor_volcano +
    geom_text_repel(data = screen$norm_tumor_top, 
                    aes(label = response, 
                        fontface = fontface, 
                        alpha = alpha), 
                    size = 2.6) + 
    geom_vline(xintercept = 0, 
               linetype = 'dashed') + 
    scale_alpha_manual(values = c(low = 0.4, 
                                  high = 1)) +
    guides(alpha = 'none') + 
    labs(subtitle = paste0('total genes: n = ', 
                           length(screen$variables), 
                           ', ', 
                           screen$norm_tumor_volcano$labels$tag), 
         tag = screen$norm_tumor_n)
  
# Correlation with the T cell content ------
  
  insert_msg('Correlation with the CD8+ T cell levels')
  
  ## Spearman correlation (CD8+ T cells fractions are not normally distributed)
  
  screen$corr_test <- screen$variables %>% 
    map(~c('T cell CD8+', .x)) %>% 
    map_dfr(~correlate_variables(screen$analysis_tbl$t_cells, 
                                 variables = .x, 
                                 what = 'correlation', 
                                 type = 'spearman', 
                                 ci = TRUE, 
                                 pub_styled = FALSE)) %>% 
    mutate(highlight = ifelse(variable2 %in% globals$cxcl_genes, 
                              paste0('<span style = "color: #000000"><b><em>', 
                                     variable2, 
                                     '</em></b></span>'), 
                              paste0('<span style = "color: #808080"><em>', 
                                     variable2, '</em></span>')))
  
  screen$corr_labs <- set_names(screen$corr_test$highlight, 
                                screen$corr_test$variable2)
  
  ## Forest plot with the correlation coefficients
  
  screen$corr_forest <- 
    plot_top(screen$corr_test, 
             regulation_variable = 'estimate', 
             label_variable = 'variable2', 
             p_variable = 'p_adjusted', 
             signif_level = 0.05, 
             regulation_level = 0, 
             lower_ci_variable = 'lower_ci', 
             upper_ci_variable = 'upper_ci', 
             top_regulated = 50, 
             plot_title = paste('Chemokine genes and CD8+ T cells, TCGA KIRC'), 
             plot_subtitle = paste('Correlation with', 
                                   'CD8+ T cell QuanTIseq estimates, n =', 
                                   nrow(screen$analysis_tbl$t_cells)), 
             x_lab = expression(rho * ', 95% CI'), 
             cust_theme = globals$common_theme) + 
    scale_y_discrete(labels = screen$corr_labs) + 
    theme(axis.text.y = element_markdown())
    
# END ------
  
  insert_tail()