# Plotting of the results of metabolic modeling

  insert_head()

# container ------

  meta_plots <- list()

# globals ------

  insert_msg('Classification of the common reactions')

  ## numbers of regulated reactions

  meta_plots$n_total <- meta$models %>%
    map(components, 'reactions') %>%
    map(length)

  meta_plots$n_reactions <- meta$models %>%
    map(components, 'regulation') %>%
    map(mutate, p_adjusted = p.adjust(p_value, 'BH')) %>% 
    map(filter, p_adjusted < 0.05) %>%
    map(mutate,
        reg_sign = ifelse(fold_reg >1, 'activated', 'inhibited'),
        reg_sign = factor(reg_sign, c('activated', 'inhibited'))) %>%
    map(count, reg_sign)

  ## labels with total and regulated reaction numbers

  meta_plots$n_labs <- meta_plots$n_reactions %>%
    map(~map2_chr(.x[[1]], .x[[2]], paste, sep = ': n = ') %>%
          paste(collapse = ', ')) %>%
    map2(meta_plots$n_total, .,
         paste, sep = ', ') %>%
    map(~paste('total reactions: n =', .x))

# Diagnostic plots: errors -------

  insert_msg('Plots of errors')

  meta_plots$errors <- meta$models %>%
    map(plot,
        type = 'errors',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohort_lexicon$label,
         ~.x +
           labs(title = .y,
                subtitle = 'Fold-regulation and regulation errors'))

# General regulation plots ------

  insert_msg('General regulation plots')

  meta_plots$regulation <- meta$models %>%
    map(plot,
        type = 'regulation',
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohort_lexicon$label,
         ~.x +
           theme(panel.grid.major.x = element_blank()) +
           labs(title = .y,
                x = 'Reactions',
                y = expression('log'[2] * ' fold-regulation, chemox high vs low'))) %>%
    map2(., meta_plots$n_labs,
         ~.x +
           labs(subtitle = .y))

# Top regulated reactions ------

  insert_msg('Top regulated reactions')

  meta_plots$top <- meta$models %>%
    map(plot,
        type = 'top',
        n_top = 10,
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohort_lexicon$label,
         ~.x +
           labs(title = .y,
                subtitle = 'Top regulated reactions'))

# Assigning the common regulated reactions to subsystemts ------
  
  insert_msg('Assigning the common reactions to subsystems')
  
  meta_plots$cmm_reactions <- sig$recon_lexicon %>% 
    select(react_id, 
           variable, 
           label) %>% 
    set_names(c('react_id', 'subsystem', 'sub_label')) %>% 
    filter(react_id %in% meta$cmm_reactions)
  
# Forest plots for the common regulated genes -------

  insert_msg('Forest plots for the common regulated genes')

  meta_plots$common <- meta$models %>%
    map(plot,
        type = 'forest',
        relevant.reactions = meta$cmm_reactions,
        cust_theme = globals$common_theme) %>%
    map2(., globals$cohort_lexicon$label,
         ~.x +
           labs(title = .y,
                subtitle = 'Common differentially regulate reactions'))

  ## adding the faceting information

  for (i in names(meta_plots$common)) {

    meta_plots$common[[i]]$data <-
      left_join(meta_plots$common[[i]]$data,
                meta_plots$cmm_reactions,
                by = 'react_id')

  }

  meta_plots$common <- meta_plots$common %>%
    map(~.x +
          facet_grid(sub_label ~ .,
                     space = 'free',
                     scales = 'free'))

# Numbers of up- and downregulated reactions per subsystem ------
  
  insert_msg('Numbers of reactions regulated per subsystem')
  
  ## subsystem lexicon
  
  meta_plots$subsystems <-  sig$recon_lexicon %>% 
    select(react_id, 
           variable, 
           label) %>% 
    set_names(c('react_id', 'subsystem', 'sub_label'))
  
  ## regulation estimates
  
  meta_plots$regulation_tbl <- meta$models %>%
    map(components, 'regulation') %>% 
    map(mutate, p_adjusted = p.adjust(p_value, 'BH')) %>% 
    map(left_join, meta_plots$subsystems, by = 'react_id')
  
  ## total reactions per subsystem
  ## up- and downregulated reactions
  
  meta_plots$subsystem_n[c('total', 
                           'upregulated', 
                           'downregulated')] <- 
    list(meta_plots$regulation_tbl,
         meta_plots$regulation_tbl %>% 
           map(filter, 
               fold_reg > 1, p_adjusted < 0.05), 
         meta_plots$regulation_tbl %>% 
           map(filter, 
               fold_reg < 1, p_adjusted < 0.05)) %>% 
    map(~map(.x, count, subsystem))
  
  meta_plots$subsystem_n <- transpose(meta_plots$subsystem_n) %>% 
    map(reduce, full_join, by = 'subsystem') %>% 
    map(set_names, 
        c('subsystem', 'n_total', 'n_upregulated', 'n_downregulated')) %>% 
    map(~map_dfc(.x, ~ifelse(is.na(.x), 0, .x)))
  
  ## calculating the fractions, adding subsystem labels
  
  meta_plots$subsystem_n <- meta_plots$subsystem_n %>% 
    map(mutate, 
        sub_label = exchange(subsystem, 
                             dict = meta_plots$subsystems, 
                             key = 'subsystem', 
                             value = 'sub_label'), 
        perc_upregulated = n_upregulated/n_total * 100, 
        perc_downregulated= n_downregulated/n_total * 100)
  
  ## calculating the subsystem regulation significance
  ## for subsystem enrichment with the Fisher exact test
  
  meta_plots$subsystem_test$up <- 
    meta_plots$subsystem_n %>% 
    map(dlply, 
        'subsystem', 
        function(x) matrix(c(x$n_upregulated[1], 
                             x$n_total[1] - x$n_upregulated[1], 
                             0, 
                             x$n_total[1]), 
                           nrow = 2)) %>% 
    map(~map(.x, 
             fisher_test, 
             alternative = 'greater', 
             detailed = TRUE) %>% 
          compress(names_to = 'subsystem') %>% 
          mutate(p_adjusted = p.adjust(p, 'BH')))
  
  meta_plots$subsystem_test$down <- 
    meta_plots$subsystem_n %>% 
    map(dlply, 
        'subsystem', 
        function(x) matrix(c(x$n_downregulated[1], 
                             x$n_total[1] - x$n_downregulated[1], 
                             0, 
                             x$n_total[1]), 
                           nrow = 2)) %>% 
    map(~map(.x, 
             fisher_test, 
             alternative = 'greater', 
             detailed = TRUE) %>% 
          compress(names_to = 'subsystem') %>% 
          mutate(p_adjusted = p.adjust(p, 'BH')))
  
  ## appending the regulated reaction count tables
  
  meta_plots$subsystem_n <- 
    map2(meta_plots$subsystem_n, 
         meta_plots$subsystem_test$up, 
         ~left_join(.x, 
                    set_names(.y[c('subsystem', 'p_adjusted')], 
                              c('subsystem', 'p_adjusted_upregulated')), 
                    by = 'subsystem'))
  
  meta_plots$subsystem_n <- 
    map2(meta_plots$subsystem_n, 
         meta_plots$subsystem_test$down, 
         ~left_join(.x, 
                    set_names(.y[c('subsystem', 'p_adjusted')], 
                              c('subsystem', 'p_adjusted_downregulated')),
                    by = 'subsystem'))
  
# Plotting the significantly regulated subsystems -------
  
  insert_msg('Plots of the regulated reaction counts')
  
  ## upregulated reactions
  
  meta_plots$subsystem_plots$upregulated <- 
    list(data = meta_plots$subsystem_n %>% 
         map(filter, p_adjusted_upregulated < 0.05), 
       title = globals$cohort_lexicon$label) %>% 
    pmap(function(data, title) data %>% 
           ggplot(aes(x = perc_upregulated, 
                      y = reorder(sub_label, perc_upregulated))) + 
           geom_bar(stat = 'identity', 
                    fill = 'firebrick3', 
                    color = 'black') + 
           labs(title = title, 
                subtitle = 'Significantly activated metabolic subsystems', 
                x = 'activated reactions, % of subsystem'))
  
  ## inhibited reactions
  
  meta_plots$subsystem_plots$downregulated <- 
    list(data = meta_plots$subsystem_n %>% 
           map(filter, p_adjusted_downregulated < 0.05), 
         title = globals$cohort_lexicon$label) %>% 
    pmap(function(data, title) data %>% 
           ggplot(aes(x = perc_downregulated, 
                      y = reorder(sub_label, perc_downregulated))) + 
           geom_bar(stat = 'identity', 
                    fill = 'steelblue3', 
                    color = 'black') + 
           labs(title = title, 
                subtitle = 'Significantly inhibited metabolic subsystems', 
                x = 'inhibited reactions, % of subsystem'))
  
  ## common adjustment
  
  meta_plots$subsystem_plots <- meta_plots$subsystem_plots %>% 
    map(~map(.x, ~.x + 
               globals$common_theme + 
               theme(axis.title.y = element_blank())))

# END -----
  
  rm(i)

  insert_tail()
