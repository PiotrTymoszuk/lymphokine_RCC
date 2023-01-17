# Paper tables and supplementary tables

  insert_head()
  
# container ------
  
  tables <- list()
  suppl_tables <- list()
  
# Table 1: cohort characteristic ------
  
  insert_msg('Table 1: cohort characteristic')
  
  tables$cohort <- cohort$ft_table %>%
    filter(!Variable %in% c('Race',
                            'Laterality',
                            'Pathological stage',
                            'Sarcomatoid differentiation',
                            'Rhabdoid differentiation',
                            'Tumor shrinkage',
                            'Clinical benefit',
                            'IMDC risk group',
                            'Neoadjuvant treatment')) %>%
    mutate(Variable = ifelse(Variable == 'OS, days',
                             'Observation time, days',
                             Variable)) %>%
    map_dfc(stri_replace_all,
            fixed = 'Complete',
            replacement = 'complete') %>%
    map_dfc(stri_replace_all,
            fixed = 'Range',
            replacement = 'range') %>%
    map_dfc(stri_replace_all,
            regex = 'Mean.*\\nMedian\\s{1}=\\s{1}',
            replacement = '') %>%
    as_mdtable(label = 'table_1_cohort_features',
               ref_name = 'cohort',
               caption = paste('Characteristic of the study cohorts.', 
                               'Numeric variables are presented as medians', 
                               'with interquartile ranges and ranges.', 
                               'Categorical variables are presented as', 
                               'percentage and total number within', 
                               'the complete observation set.'))
  
# Supplementary Table S1: listing of significantly regulated Reactome signs ------
  
  insert_msg('Table S1: Reactome signatures')  
  
  suppl_tables$reactome <- clust_reactome$significant %>% 
    map(arrange, - estimate) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = exchange(cohort, 
                             dict = globals$cohort_lexicon, 
                             key = 'cohort', 
                             value = 'label')) %>% 
    select(cohort, sign_label,
           estimate, lower_ci, upper_ci, p_adjusted) %>% 
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>%
    set_names(c('Cohort',
                'Reactome signature',
                'delta ssGSEA score',
                'lower CI',
                'upper CI',
                'pFDR')) %>% 
    as_mdtable(label = 'table_s1_reactome',
               ref_name = 'reactome',
               caption = paste('Reactome signatures differentially regulated', 
                               'between the lymphokine clusters', 
                               '(chemox high vs low). The table is', 
                               'available as a supplementary Excel file.'))

# Supplementary Table S2: numbers of differentially regulated genes -----
  
  insert_msg('Table S2: numbers of differentially regulated genes')
  
  suppl_tables$n_dge <- clust_dge$significant %>% 
    map(count, regulation) %>% 
    map(~tibble(`N downregulated` = .x$n[2], 
                `N upregulated` = .x$n[1])) %>% 
    compress(names_to = 'Cohort') %>% 
    mutate(Cohort = exchange(Cohort, 
                             dict = globals$cohort_lexicon, 
                             key = 'cohort', 
                             value = 'label')) %>% 
    select(Cohort, `N downregulated`, `N upregulated`) %>%
    as_mdtable(label = 'table_s2_number_dge',
               ref_name = 'n_dge',
               caption = paste('Numbers of genes differentially regulated', 
                               'between the lymphokine clusters', 
                               '(chemox high vs low).'))
  
# Supplementary Table S3: listing of differentially regulated genes -------
  
  insert_msg('Table S3: differentially regulated gene list')
  
  suppl_tables$dge_genes <- clust_dge$significant %>%
    map(mutate, 
        entrez_id = as.character(entrez_id)) %>%
    map(arrange, - estimate) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = exchange(cohort, 
                             dict = globals$cohort_lexicon, 
                             key = 'cohort', 
                             value = 'label')) %>% 
    select(cohort, regulation,
           gene_symbol, entrez_id,
           estimate, lower_ci, upper_ci, p_adjusted) %>%
    map_dfc(function(x) if(is.numeric(x)) signif(x, 3) else x) %>%
    set_names(c('Cohort',
                'Regulation',
                'Gene symbol',
                'Entrez ID',
                'Log2 fold-regulation estimate',
                'lower CI',
                'upper CI',
                'pFDR')) %>%
    as_mdtable(label = 'table_s3_dge',
               ref_name = 'dge_gene',
               caption = paste('Genes differentially regulated between', 
                               'the lymphokine clusters (chemox high vs low).', 
                               'The table is', 
                               'available as a supplementary Excel file.'))
  
# Supplementary Table S4: signaling pathways -------
  
  insert_msg('Table S4: signaling pathways')
  
  suppl_tables$spia <- dge_spia$test_results %>%
    map(filter, pGFdr < 0.05) %>%
    map(arrange, -tA) %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = exchange(cohort,
                             dict = globals$cohort_lexicon, 
                             key = 'cohort', 
                             value = 'label')) %>%
    select(cohort,
           ID,
           Name,
           tA,
           pGFdr) %>%
    set_names(c('Cohort',
                'KEGG pathway ID',
                'KEGG pathway name',
                'Magnitude of regulation, tA',
                'pFDR')) %>%
    as_mdtable(label = 'table_s5_spia',
               ref_name = 'spia',
               caption = paste('KEGG (Kyoto Encyclopedia of Genes and', 
                               'Genomes) signaling pathways predicted to be', 
                               'significantly regulated between the', 
                               'lymphokine clusters (chemox high vs low)', 
                               'by the SPIA tool.', 
                               'The table is available as a supplementary', 
                               'Excel file.'))
  
# Supplementary Table S5: metabolic pathways ------
  
  insert_msg('Table S5: predicted regulation of metabolic reactions')
  
  suppl_tables$reacts <- meta$models %>%
    map(components, 'regulation') %>%
    map(filter, p_value < 0.05) %>%
    map(mutate,
        react_name = annotate_bigg(react_id,
                                   annotation_db = biggrExtra::Recon2D),
        react_name = unlist(react_name)) %>%
    map(left_join, 
        sig$recon_lexicon[c('react_id', 'subsystem')], 
        by = 'react_id') %>% 
    map(arrange, subsystem, -fold_reg) %>% 
    compress(names_to = 'cohort') %>% 
    mutate(cohort = exchange(cohort, 
                             dict = globals$cohort_lexicon,
                             key = 'cohort', 
                             value = 'label')) %>%
    select(cohort,
           subsystem, 
           react_id,
           react_name,
           fold_reg,
           lower_ci,
           upper_ci,
           p_value) %>%
    set_names(c('Cohort',
                'Subsystem', 
                'BiGG reaction ID',
                'BiGG reaction name',
                'Activity fold regulation',
                'lower CI',
                'upper CI',
                'pFDR')) %>%
    as_mdtable(label = 'table_s6_reactions',
               ref_name = 'reacts',
               caption = paste('Recon model metabolic reactions predicted to be', 
                               'significantly modulated between the', 
                               'lymphokine clusters (chemox high vs low).', 
                               'The table is available as a', 
                               'supplementary Excel file.'))
  
  
# Saving the Supplementary tables -----
  
  insert_msg('Saving the supplementary tables')
  
  ## cover page with the legend captions
  
  suppl_tables$cover <- suppl_tables %>% 
    map_chr(attr, 'caption') %>% 
    tibble(`Supplementary Table` = paste0('Table S', 1:length(suppl_tables)), 
           Caption = .)
  
  ## saving the supplementary tables
  
  suppl_tables[c('cover', names(suppl_tables)[names(suppl_tables) != 'cover'])] %>% 
    set_names(c('Cover', paste0('Table S', 1:(length(suppl_tables) - 1)))) %>%
    write_xlsx('./paper/supplementary_tables.xlsx')
  
  # END ------
  
  insert_tail()
  