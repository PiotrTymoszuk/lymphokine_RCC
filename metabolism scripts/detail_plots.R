# Detailed plots for OxPhos, TCA, FAO and TRP metabolism

  insert_head()
  
# container -------
  
  meta_details <- list()
  
# analysis globals -------
  
  insert_msg('Analysis globals')
  
  ## subsystems of interest
  
  meta_details$subsystems <- sig$recon_lexicon %>% 
    select(react_id, 
           variable, 
           label) %>% 
    set_names(c('react_id', 'subsystem', 'sub_label')) %>% 
    filter(subsystem %in% c('oxidative_phosphorylation', 
                            'citric_acid_cycle', 
                            'fatty_acid_oxidation', 
                            'tryptophan_metabolism', 
                            'glycolysis_gluconeogenesis'))
  
  ## regulation estimates for the sub-systems of interest
  
  meta_details$analysis_tbl <- meta$models %>%
    map(components, 'regulation') %>%
    map(mutate, 
        log_fold_reg = log2(fold_reg), 
        p_adjusted = p.adjust(p_value, 'BH'), 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(fold_reg > 1, 
                                   'activated', 'inhibited')), 
        regulation = factor(regulation, 
                            c('activated', 'inhibited', 'ns'))) %>% 
    map(~left_join(meta_details$subsystems, 
                   .x, 
                   by = 'react_id'))
  
# Sign plots for the fatty acid oxidation -------
  
  insert_msg('Fatty acid oxidation')
  
  ## numbers of up- and down-regulated reactions
  
  meta_details$fao_n <- meta_details$analysis_tbl %>% 
    map(filter, subsystem == 'fatty_acid_oxidation') %>% 
    map(count, regulation, .drop = FALSE) %>% 
    map(mutate, 
        n_total = sum(n)) %>% 
    map(~paste0('total: n = ', 
                .x[[3]][1], 
                ', activated: n = ', 
                .x[[2]][1], 
                ', inhibited: n = ', 
                .x[[2]][2]))
  
  ## sign plots 
  
  meta_details$fao_sign_plots <- 
    list(data = meta_details$analysis_tbl %>% 
           map(filter, subsystem == 'fatty_acid_oxidation'), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle =  meta_details$fao_n) %>%
    pmap(plot_sign, 
         regulation_variable = 'log_fold_reg', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         cust_theme = globals$common_theme, 
         y_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         x_lab = 'Fatty acid oxidation reaction', 
         show_trend = FALSE) %>% 
    map(~.x + 
          scale_fill_manual(values = c(upregulated = 'firebrick', 
                                       downregulated = 'steelblue', 
                                       ns = 'gray60'), 
                            labels = c(upregulated = 'activated', 
                                       downregulated = 'inhibited', 
                                       ns = 'ns'), 
                            name = 'Regulation\nchemox high vs low') + 
          theme(plot.tag = element_blank()))
  
# Sign plots for the tryptophan metabolism -------
  
  insert_msg('Tryptophan metabolism')
  
  ## numbers of up- and down-regulated reactions
  
  meta_details$trp_n <- meta_details$analysis_tbl %>% 
    map(filter, subsystem == 'tryptophan_metabolism') %>% 
    map(count, regulation, .drop = FALSE) %>% 
    map(mutate, 
        n_total = sum(n)) %>% 
    map(~paste0('total: n = ', 
                .x[[3]][1], 
                ', activated: n = ', 
                .x[[2]][1], 
                ', inhibited: n = ', 
                .x[[2]][2]))
  
  ## sign plots 
  
  meta_details$trp_sign_plots <- 
    list(data = meta_details$analysis_tbl %>% 
           map(filter, subsystem == 'tryptophan_metabolism'), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle =  meta_details$trp_n) %>%
    pmap(plot_sign, 
         regulation_variable = 'log_fold_reg', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         cust_theme = globals$common_theme, 
         y_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         x_lab = 'Tryptophan metabolism reaction', 
         show_trend = FALSE) %>% 
    map(~.x + 
          scale_fill_manual(values = c(upregulated = 'firebrick', 
                                       downregulated = 'steelblue', 
                                       ns = 'gray60'), 
                            labels = c(upregulated = 'activated', 
                                       downregulated = 'inhibited', 
                                       ns = 'ns'), 
                            name = 'Regulation\nchemox high vs low') + 
          theme(plot.tag = element_blank()))
  
# Detailed plots for the IDO - TRYCATS pathway -------
  
  insert_msg('Detailed plots for the IDO/TRYCATS pathway')
  
  ## reactions, including BH4 biosynthesis pathway
  
  meta_details$trycats_reactions <- c('R_5HTRPDOX', 
                                      'R_TRPO2', 
                                      'R_KYN3OX',
                                      'R_KYNAKGAT', 
                                      'R_KYN', 
                                      'R_3HAO')
  
  ## regulation estimates
  
  meta_details$trycats_tbl <- meta$models %>%
    map(components, 'regulation') %>%
    map(mutate, 
        log_fold_reg = log2(fold_reg), 
        p_adjusted = p.adjust(p_value, 'BH'), 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(fold_reg > 1, 
                                   'activated', 'inhibited')), 
        regulation = factor(regulation, 
                            c('activated', 'inhibited', 'ns'))) %>% 
    map(filter, 
        react_id %in% meta_details$trycats_reactions) %>% 
    map(mutate, 
        react_name = annotate_bigg(react_id), 
        react_name = stri_replace_all(react_name, 
                                      regex = '\\(.*\\)', 
                                      replacement = ''), 
        react_id = factor(react_id, 
                          meta_details$trycats_reactions), 
        react_id = droplevels(react_id))
  
  ## regulation plots
  
  meta_details$trycats_detail_plots <- 
    list(data = meta_details$trycats_tbl, 
         title = globals$cohort_lexicon$label) %>% 
    pmap(function(data, title) data %>% 
           ggplot(aes(x = log_fold_reg, 
                      y = reorder(react_name, 
                                  -as.numeric(react_id)), 
                      fill = regulation)) + 
           geom_bar(stat = 'identity', 
                    color = 'black') + 
           scale_fill_manual(values = c(activated = 'firebrick', 
                                        inhibited = 'steelblue', 
                                        ns = 'gray60'), 
                             name = 'Regulation\nchemox high vs chemox low') + 
           globals$common_theme + 
           theme(axis.title.y = element_blank()) + 
           labs(title = title, 
                subtitle = 'IDO/TRACATS pathway reactions', 
                x = expression('log'[2] * ' fold regulation, chemox high vs low')))
  
# Sign plots for the glycolysis/gluconeogenesis subsystem -------
  
  insert_msg('Glycolysis/gluconeogenesis')
  
  ## numbers of up- and down-regulated reactions
  
  meta_details$glyco_n <- meta_details$analysis_tbl %>% 
    map(filter, subsystem == 'glycolysis_gluconeogenesis') %>% 
    map(count, regulation, .drop = FALSE) %>% 
    map(mutate, 
        n_total = sum(n)) %>% 
    map(~paste0('total: n = ', 
                .x[[3]][1], 
                ', activated: n = ', 
                .x[[2]][1], 
                ', inhibited: n = ', 
                .x[[2]][2]))
  
  ## sign plots 
  
  meta_details$glyco_sign_plots <- 
    list(data = meta_details$analysis_tbl %>% 
           map(filter, subsystem == 'glycolysis_gluconeogenesis'), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle =  meta_details$glyco_n) %>%
    pmap(plot_sign, 
         regulation_variable = 'log_fold_reg', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         cust_theme = globals$common_theme, 
         y_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         x_lab = 'Glycolysis/gluconeogenesis reaction', 
         show_trend = FALSE) %>% 
    map(~.x + 
          scale_fill_manual(values = c(upregulated = 'firebrick', 
                                       downregulated = 'steelblue', 
                                       ns = 'gray60'), 
                            labels = c(upregulated = 'activated', 
                                       downregulated = 'inhibited', 
                                       ns = 'ns'), 
                            name = 'Regulation\nchemox high vs low') + 
          theme(plot.tag = element_blank()))
  
# oxidative phosphorylation -----
  
  insert_msg('Ovidative phosphorylation')
  
  ## plotting tables, selecting the mitochondrial OxPhos complex reactions
  
  meta_details$oxphos_tbl <- meta$models %>%
    map(components, 'regulation') %>% 
    map(mutate, 
        log_fold_reg = log2(fold_reg), 
        p_adjusted = p.adjust(p_value, 'BH'), 
        regulation = ifelse(p_adjusted >= 0.05, 
                            'ns', 
                            ifelse(fold_reg > 1, 
                                   'activated', 'inhibited'))) %>% 
    map(filter, 
        react_id %in% c('R_ATPS4m',
                        'R_CYOOm2',
                        #'R_CYOOm3', ## superoxide production
                        'R_CYOR_u10m',
                        'R_FADH2ETC',
                        'R_NADH2_u10m')) %>% 
    map(mutate, 
        complex = car::recode(react_id,
                              "'R_NADH2_u10m' = 'I';
                              'R_FADH2ETC' = 'II';
                              'R_CYOR_u10m' = 'III';
                              'R_CYOOm2' = 'IV';
                              'R_CYOOm3' = 'IV';
                              'R_ATPS4m' = 'V'"),
        complex = factor(complex, c('I', 'II', 'III', 'IV', 'V')), 
        regulation = ifelse(is.na(regulation), 'ns', regulation),
        regulation = factor(regulation, 
                            c('activated', 'inhibited', 'ns')))
  
  ## bar plots with the regulation signs
  
  meta_details$oxphos_sign_plots <- 
    list(data = meta_details$oxphos_tbl, 
         title = globals$cohort_lexicon$label) %>% 
    pmap(function(data, title) data %>% 
           ggplot(aes(x = log_fold_reg, 
                      y = complex, 
                      fill = regulation)) + 
           geom_bar(stat = 'identity', 
                    color = 'black') + 
           scale_fill_manual(values = c(activated = 'firebrick', 
                                        inhibited = 'steelblue', 
                                        ns = 'gray60'), 
                             name = 'Regulation\nchemox high vs low') + 
           globals$common_theme + 
           labs(title = title, 
                subtitle = 'Activity of OxPhos complexed, chemox high vs low', 
                x = expression('log'[2] * ' fold regulation, chemox high vs low'), 
                y = 'OxPhos complex'))
  
# Sign plots for the citric acid cycle subsystem ---------
  
  insert_msg('Citric acid cycle')
  
  ## numbers of up- and down-regulated reactions
  
  meta_details$tca_n <- meta_details$analysis_tbl %>% 
    map(filter, subsystem == 'citric_acid_cycle') %>% 
    map(count, regulation, .drop = FALSE) %>% 
    map(mutate, 
        n_total = sum(n)) %>% 
    map(~paste0('total: n = ', 
                .x[[3]][1], 
                ', activated: n = ', 
                .x[[2]][1], 
                ', inhibited: n = ', 
                .x[[2]][2]))
  
  ## sign plots 
  
  meta_details$tca_sign_plots <- 
    list(data = meta_details$analysis_tbl %>% 
           map(filter, subsystem == 'citric_acid_cycle'), 
         plot_title = globals$cohort_lexicon$label, 
         plot_subtitle =  meta_details$tca_n) %>%
    pmap(plot_sign, 
         regulation_variable = 'log_fold_reg', 
         p_variable = 'p_adjusted', 
         signif_level = 0.05, 
         regulation_level = 0, 
         cust_theme = globals$common_theme, 
         y_lab = expression('log'[2] * ' fold regulation, chemox high vs low'), 
         x_lab = 'Citric acid cycle reaction', 
         show_trend = FALSE) %>% 
    map(~.x + 
          scale_fill_manual(values = c(upregulated = 'firebrick', 
                                       downregulated = 'steelblue', 
                                       ns = 'gray60'), 
                            labels = c(upregulated = 'activated', 
                                       downregulated = 'inhibited', 
                                       ns = 'ns'), 
                            name = 'Regulation\nchemox high vs low') + 
          theme(plot.tag = element_blank()))
  
# Detailed plots for the citric acid cycle reactions -------
  
  insert_msg('Detailed plots for the TCA reactions')
  
  ## plotting data
  
  meta_details$tca_plot_tbl <- meta_details$analysis_tbl %>% 
    map(filter, 
        subsystem == 'citric_acid_cycle', 
        react_id != 'R_ICDHy') %>% ## duplicated reaction
    map(mutate, 
        react_name = annotate_bigg(react_id), 
        react_name = ifelse(is.na(react_name), 
                            car::recode(react_id, 
                                        "'R_ICDHyrm' = 'Isocitrate dehydrogenase NADP, mitochondrial'; 
                                        'R_ACONTm' = 'Aconitase hydratase, mitochondrial'; 
                                        'R_r0081' = 'Oxalosuccinate:NADP+ oxidoreductase'"), 
                            react_name), 
        react_name = ifelse(react_name == 'R_r0081', 
                            'Oxalosuccinate:NADP+ oxidoreductase', 
                            react_name), 
        react_name = ifelse(react_id == 'R_ACONTm', 
                            'Aconitase hydratase, mitochondrial', 
                            react_name), 
        react_id = factor(react_id, 
                          c('R_ACONT', 
                            'R_ACONTm', 
                            'R_ICDHxm', 
                            'R_ICDHy', 
                            'R_ICDHyp', 
                            'R_ICDHyrm', 
                            'R_r0081', 
                            'R_AKGDm', 
                            'R_SUCOAS1m', 
                            'R_SUCOASm', 
                            'R_SUCD1m', 
                            'R_FUM', 
                            'R_FUMm', 
                            'R_MDH', 
                            'R_MDHm', 
                            'R_CSm'))) %>% 
    map(filter, !is.na(react_id))
  
  ## plots 
  
  meta_details$tca_detail_plots <- 
    list(data = meta_details$tca_plot_tbl, 
         title = globals$cohort_lexicon$label) %>% 
    pmap(function(data, title) data %>% 
           ggplot(aes(x = log_fold_reg, 
                      y = reorder(react_name, 
                                  -as.numeric(react_id)), 
                      fill = regulation)) + 
           geom_bar(stat = 'identity', 
                    color = 'black') + 
           scale_fill_manual(values = c(activated = 'firebrick', 
                                        inhibited = 'steelblue', 
                                        ns = 'gray60'), 
                             name = 'Regulation\nchemox high vs chemox low') + 
           globals$common_theme + 
           theme(axis.title.y = element_blank()) + 
           labs(title = title, 
                subtitle = 'Citric acid reactions', 
                x = expression('log'[2] * ' fold regulation, chemox high vs low')))
  
# END -------
  
  insert_tail()