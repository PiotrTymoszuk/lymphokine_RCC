# This script provides tools for project data analyses

# tools ----

  library(plyr)
  library(tidyverse)
  library(rlang)
  library(trafo)
  library(ggrepel)
  library(survival)
  library(coxExtensions)
  library(igraph)
  library(network)
  library(ggnetwork)

# annotation functions -----

  translate_gene <- function(id,
                             id_type = 'gene_symbol',
                             output_type = 'entrez_id',
                             dictionary = tcga$annotation) {

    ## gets gene identifier of interest

    naming_vec <- dictionary[[output_type]] %>%
      set_names(dictionary[[id_type]])

    naming_vec[id]

  }

  mm_inch <- function(x) x * 0.0393700787

# bubble plots of correlation with immune infiltrates ------
  
  plot_immune_bubble <- function(correlation_data, 
                                 plot_title = NULL, 
                                 variable_lexicon = globals$quantiseq_lexicon) {
    
    ## bubble plot of correlation coefficients
    ## significant effects are labeled in bold
    
    ggplot(correlation_data,
           aes(x = variable1,
               y = variable2,
               size = abs(estimate),
               fill = estimate)) +
      geom_point(shape = 21) +
      geom_text(aes(label = plot_lab,
                    fontface = significant, 
                    alpha = significant),
                size = 2.75,
                hjust = 0.5,
                vjust = -1.5) +
      scale_y_discrete(labels = set_names(variable_lexicon[['label']], 
                                          variable_lexicon[['variable']])) +
      scale_fill_gradient2(low = 'steelblue',
                           mid = 'white',
                           high = 'firebrick',
                           midpoint = 0,
                           limits = c(-1, 1),
                           name = expression(rho)) +
      scale_size_area(limits = c(0, 1),
                      max_size = 6,
                      name = expression('abs('*rho*')')) +
      scale_x_discrete(limits = globals$cxcl_genes) +
      scale_alpha_manual(values = c(1, 0.4)) + 
      globals$common_theme +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(face = 'italic')) +
      labs(title = plot_title,
           subtitle = paste('n =', correlation_data$n[[1]]))
    
  }
  
  plot_simplified_bubble <- function(correlation_data, 
                                     plot_title = NULL, 
                                     max_size = 3.5, 
                                     variable_lexicon = globals$quantiseq_lexicon) {
    
    correlation_data %>% 
      ggplot(aes(x = variable1, 
                 y = variable2, 
                 fill = estimate_recoded, 
                 size = abs(estimate_recoded))) + 
      geom_point(shape = 21) + 
      scale_size_area(max_size = max_size, 
                      limits = c(0, 1), 
                      name = expression('abs(' * rho * ')')) + 
      scale_x_discrete(limits = globals$cxcl_genes) + 
      scale_y_discrete(labels = set_names(variable_lexicon[['label']], 
                                          variable_lexicon[['variable']]), 
                       limits = rev(variable_lexicon[['variable']])) +
      scale_fill_gradient2(low = 'steelblue',
                           mid = 'white',
                           high = 'firebrick',
                           midpoint = 0,
                           limits = c(-1, 1),
                           name = expression(rho)) + 
      globals$common_theme + 
      theme(axis.title = element_blank(), 
            axis.text.x = element_text(face = 'italic')) + 
      labs(title = plot_title, 
           subtitle = paste('n =', correlation_data[['n']][1]))
    
  }
  
# UMAP and MDS ------
  
  plot_umap <- function(data, 
                        overlay = 'CXCL9', 
                        plot_title = overlay, 
                        plot_subtitle = NULL, 
                        plot_tag = NULL, 
                        fill_title = 'Z-score', 
                        center = 'mean', 
                        point_size = 2, 
                        point_alpha = 1, 
                        point_circle = 'none', 
                        low_color = 'steelblue3', 
                        mid_color = 'gray80', 
                        high_color = 'firebrick3', 
                        mid_point = 0, ...) {
    
    ## plots UMAP dimensions with the overlay variable
    ## coded by point color. If center is 'mean' or 'median', Z-scores
    ## with the respective centering are calculated
    
    center <- match.arg(center[1], c('none', 'mean', 'median'))
    
    if(center != 'none') {
      
      center_fun <- switch(center, 
                           mean = function(x) mean(x, na.rm = TRUE), 
                           median = function(x) median(x, na.rm = TRUE))
      
      data <- data %>% 
        mutate(!!overlay := scale(.data[[overlay]], 
                                  center = center_fun(.data[[overlay]]))[, 1])
      
    }
    
    data <- data[c('comp_1', 'comp_2', overlay)]
    
    if(is.null(plot_tag)) {
      
      plot_tag <- paste('\nn =', nrow(data))
      
    }
    
    ## plotting

    if(point_circle == 'none') {
      
      umap_plot <- data %>% 
        ggplot(aes(x = comp_1, 
                   y = comp_2, 
                   color = .data[[overlay]])) + 
        geom_point(shape = 16, 
                   alpha = point_alpha, 
                   size = point_size) + 
        scale_color_gradient2(low = low_color, 
                              mid = mid_color, 
                              high = high_color, 
                              midpoint = mid_point, 
                              name = fill_title, ...)
      
    } else {
      
      umap_plot <- data %>% 
        ggplot(aes(x = comp_1, 
                   y = comp_2, 
                   fill = .data[[overlay]])) + 
        geom_point(shape = 21, 
                   size = point_size, 
                   color = point_circle) + 
        scale_fill_gradient2(low = low_color, 
                             mid = mid_color, 
                             high = high_color, 
                             midpoint = mid_point, 
                             name = fill_title, ...)
      
      
    }
 
    umap_plot +  
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = 'UMAP 1', 
           y = 'UMAP 2')
    
  }
  
  plot_mds <- function(mds_object, 
                       knn_data, 
                       variable_lexicon = globals$xcell_lexicon, 
                       plot_title = NULL, 
                       plot_subtitle = NULL, 
                       plot_tag = NULL, 
                       point_size = 2, 
                       point_circle = 'none', 
                       txt_size = 2.6, 
                       nearest_only = TRUE) {
    
    ## plots the MDS layout with the genes of interest labeled
    ## along with their nearest immune neighbors
    
    ## NN labels
    
    knn_labels <- knn_data %>% 
      mutate(knn = map(knn, 
                       exchange, 
                       dict = variable_lexicon, 
                       key = 'variable', 
                       value = 'label'))
    
    knn_labels <- knn_labels$knn %>% 
      reduce(union)
    
    knn_names <- knn_data$knn %>% 
      reduce(union)
    
    ## plot data
    
    plot_tbl <- mds_object$component_tbl %>% 
      mutate(interest = ifelse(observation %in% globals$cxcl_genes, 
                               'yes', 'no'), 
             interest = factor(interest, c('no', 'yes')), 
             interest_name = ifelse(interest == 'yes', observation, NA), 
             nearest = ifelse(observation %in% knn_names, 
                              exchange(observation, 
                                       dict = variable_lexicon, 
                                       key = 'variable', 
                                       value = 'label'), 
                              NA), 
             interest_name = ifelse(is.na(interest_name), 
                                    nearest, 
                                    interest_name), 
             fontface = ifelse(interest == 'yes', 
                               'italic', 'plain'), 
             point_alpha = ifelse(is.na(interest_name), 'low', 'high'))
    
    if(nearest_only) {
      
      plot_tbl <- plot_tbl %>% 
        filter(!is.na(interest_name))
      
    }
    
    if(is.null(plot_tag)) {
      
      plot_tag <- paste('\nn =', nobs(mds_object)$variables)
      
    }
    
    ## plotting
    
    if(point_circle == 'none') {
      
      mds_plot <- plot_tbl %>% 
        ggplot(aes(x = comp_1, 
                   y = comp_2, 
                   color = interest, 
                   alpha = point_alpha)) + 
        geom_point(shape = 16, 
                   size = point_size)
      
    } else {
      
      mds_plot <- plot_tbl %>% 
        ggplot(aes(x = comp_1, 
                   y = comp_2, 
                   fill = interest, 
                   alpha = point_alpha)) + 
        geom_point(shape = 21, 
                   color = point_circle, 
                   size = point_size)
      
    }
    
    mds_plot <- mds_plot + 
      geom_text_repel(aes(label = interest_name, 
                          color = interest, 
                          fontface = fontface), 
                      size = txt_size) + 
      scale_fill_manual(values = c(no = 'steelblue4', 
                                   yes = 'firebrick4')) +
      scale_color_manual(values = c(no = 'steelblue4', 
                                    yes = 'firebrick4')) + 
      scale_alpha_manual(values = c(high = 1, 
                                    low = 0.25)) + 
      guides(fill = 'none', 
             color = 'none', 
             alpha = 'none') + 
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = 'MDS 1', 
           y = 'MDS 2')
    
    return(mds_plot)
    
  }
  
  forced_netgraph <- function(data, 
                              corr_limit = 0.4, 
                              cell_subset = globals$xcell_lexicon$variable, 
                              gene_subset = globals$cxcl_genes, 
                              weighting_order = 1, 
                              gene_color = 'coral3', 
                              cell_color = 'steelblue', 
                              txt_size = 2.75, 
                              plot_title = NULL, 
                              plot_subtitle = NULL, 
                              cust_theme = theme_blank()) {
    
    ## creates a force-directed graph of a min/max scaled Spearman's 
    ## correlation matrix
    ## for genes and cell types quantified by xCell
    ##
    ## data is a data frame with the expression levels and predicted 
    ## infiltration quantities
    ##
    ## corr_limit sets correlations below the cutoff to 0
    ##
    ## cell_subset tells which cell types will be plotted
    ## gene_subset tells which genes will be plotted
    ## 
    ## weighting_order allows for quadratic, cubic etc. transformation 
    ## of the scaled Spearman's correlation coefficient - this one codes 
    ## for node alpha
    ##
    ## gene_color and cell_color define colors of the nodes representing 
    ## genes and cells, respectively
    
    ## correlation matrix: min/max scaled Spearman's coefficient
    ## squishing the correlations
    
    corr_mtx <- cor(data, method = 'spearman')
    
    corr_mtx <- (corr_mtx - min(corr_mtx))/(max(corr_mtx) - min(corr_mtx))
    
    corr_mtx[corr_mtx < corr_limit] <- 0
    
    ## variable lexicon
    
    var_lexicon <- rbind(globals$xcell_lexicon, 
                         tibble(variable = globals$cxcl_genes, 
                                label = globals$cxcl_genes))
    
    
    
    ## setting the cell populations and genes to be plotted
    ## weighted graph with `graph_from_adjacency_matrix()`, package igraph
    
    plot_mtx <- corr_mtx[c(cell_subset, gene_subset), 
                         c(cell_subset, gene_subset)]
    
    plot_graph <- graph_from_adjacency_matrix(plot_mtx, 
                                              mode = 'undirected', 
                                              diag = FALSE, 
                                              weighted = TRUE) %>% 
      set_vertex_attr(name = 'node_type', 
                      value = c(rep('cell', length(cell_subset)), 
                                rep('gene', length(gene_subset)))) %>% 
      set_vertex_attr(name = 'font_type', 
                      value = c(rep('plain', length(cell_subset)), 
                                rep('italic', length(gene_subset)))) %>% 
      set_vertex_attr(name = 'node_lab', 
                      value = exchange(c(cell_subset, gene_subset), 
                                       dict = var_lexicon, 
                                       key = 'variable', 
                                       value = 'label'))
    
    ## plotting 
    
    if(weighting_order == 1) {
      
      alpha_title <- 'scaled \u03C1'
      
    } else {
      
      alpha_title <- paste0('scaled \u03C1 <sup>', weighting_order, '</sup>')
      
    }
    
    ggplot(plot_graph, 
           aes(x = x, y = y, xend = xend, yend = yend)) + 
      geom_edges(aes(alpha = weight^weighting_order)) + 
      geom_nodes(aes(color = node_type), 
                 size = 3) + 
      geom_nodetext_repel(aes(label = node_lab, 
                              color = node_type, 
                              fontface = font_type), 
                          size = txt_size) + 
      scale_color_manual(values = c(gene = gene_color, 
                                    cell = cell_color)) + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           color = 'feature type', 
           alpha = alpha_title) + 
      cust_theme + 
      theme(legend.title = element_markdown(size = 8), 
            legend.text = globals$common_text, 
            plot.subtitle = globals$common_text, 
            plot.title = element_text(size = 8, face = 'bold'))
    
  }
  
# Splined Cox models -----
  
  plot_spline_cox <- function(coxex_model, 
                              plot_title = NULL, 
                              color = 'steelblue', 
                              alpha = 0.25, 
                              x_lab = 'Z score expression', 
                              y_lab = 'log HR, 95% CI') {
    
    ## plots predictions of log OR for the non-linear term
    
    ## inference
    
    coefs <- summary(coxex_model$model)$coefficients %>% 
      as.data.frame %>% 
      mutate(significance = ifelse(p < 0.05, 
                                   paste('p =', signif(p, 2)), 
                                   paste0('ns (p = ', signif(p, 2), ')')), 
             term = c('linear', 'spline'), 
             stat_lab = paste0('\u03C7\u00B2(', signif(DF, 2), 
                               ') = ', signif(Chisq, 2)), 
             plot_cap = paste0(term, ': ', stat_lab, ', ', significance))
    
    plot_subtitle <- coefs$plot_cap %>% 
      paste(collapse = '\n')
    
    ## plotting data
    
    plot_tbl <- termplot(coxex_model$model, 
                         se = TRUE, 
                         plot = FALSE)[[1]] %>% 
      transmute(expression = x, 
                estimate = y, 
                se = se, 
                lower_ci = estimate + qnorm(0.025) * se, 
                upper_ci = estimate + qnorm(0.975) * se) %>% 
      as_tibble
    
    plot_tag <- paste0('\ntotal: n = ', 
                       coxex_model$model$n, 
                       ', events: n = ', 
                       coxex_model$model$nevent)
    
    ## plotting
    
    spline_plot <- plot_tbl %>% 
      ggplot(aes(x = expression, 
                 y = estimate)) +
      geom_hline(yintercept = 0, 
                 linetype = 'dashed') + 
      geom_ribbon(aes(ymin = lower_ci, 
                      ymax = upper_ci), 
                  fill = color, 
                  color = NA, 
                  alpha = alpha) + 
      geom_path(color = color) +
      globals$common_theme + 
      labs(title = plot_title, 
           subtitle = plot_subtitle, 
           tag = plot_tag, 
           x = x_lab, 
           y = y_lab)
    
    return(spline_plot)
    
  }
  
# labellers ------
  
  reactome_labeller <- function(x) {
    
    x <- tolower(x)
    
    replacements <- 
      list(c('pd 1', 'PD1'), 
           c('generation of ', ''), 
           c('flt3', 'FLT3'), 
           c('src', 'SRC'), 
           c('family ', ''), 
           c('interleukin ', 'IL'), 
           c('interferon ', 'IFN '), 
           c('b cell receptor ', ''), 
           c('cd28', 'CD28'), 
           c('nf kb', 'NF-kB'), 
           c('ca 2', 'Ca2+'), 
           c('dap12', 'DAP12'), 
           c('il10', 'IL10'), 
           c('fcgamma receptor ', ''), 
           c('fc epsilon receptor ', ''), 
           c(' favouring ', ', '), 
           c('fcgr', 'FCGR'), 
           c('bcr', 'BCR'), 
           c('nef', 'NEF'), 
           c('fceri', 'FCERI'), 
           c('fcgr3a', 'FCGR3A'),
           c('gpvi', 'GPVI'), 
           c('immunoregulatory interactions between a lymphoid and a non lymphoid cell', 
             'immunoregulation: lymphoid \u2194 non-lymphoid cell'), 
           c('leishmania parasite', 
             'leishmania'))
    
    for(i in replacements) {
      
      x <- stri_replace_all(x, fixed = i[1], replacement = i[2])
      
    }
    
    return(x)
    
  }
  
  metabolism_labeller <- function(x) {
    
    replacements <- 
      list(c('Androgen and estrogen', 'AR/E'), 
             c('Vitamin', 'Vit.'), 
             c('Glycine, serine, alanine and threonine', 
               'Gly/Ser/Ala/Thr'), 
           c('Phosphatidylinositol phosphate', 'PIP'), 
           c('Tetrahydrobiopterin', 'BH4'), 
           c('Valine, leucine, and isoleucine', 
             'Val/Leu/Ile'), 
           c('Methionine and cysteine', 'Met/Cys'), 
           c('Arginine and Proline', 'Arg/Pro'), 
           c('Tryptophan', 'TRP'), 
           c('Fructose and mannose', 'Fruc/Man'), 
           c('dibasic acid', 'DBA'), 
           c('Glycerophospholipid', 'GP'))
    
    
    for(i in replacements) {
      
      x <- stri_replace_all(x, fixed = i[1], replacement = i[2])
      
    }
    
    ifelse(x %in% c('Oxidative phosphorylation', 
                    'Citric acid cycle', 
                    'Fatty acid oxidation', 
                    'TRP metabolism'), 
           paste0('<span style = "color: #000000"><b>', 
                  x, 
                  '</b></span>'), 
           paste0('<span style = "color: #808080">', 
                  x, '</span>'))
    
  }
  
  tca_labeller <- function(x) {
    
    replacements <- 
      list(c('Aconitate hydratase', 'ACON'), 
           c('Aconitase hydratase, mitochondrial', 'ACON mito'), 
           c('Isocitrate dehydrogenase  NAD', 'IDH/NAD'), 
           c('Isocitrate dehydrogenase  NADP', 'IDH/NADP'), 
           c('Isocitrate dehydrogenase NADP, mitochondrial', 'IDH/NADP mito'), 
           c('Oxalosuccinate:NADP+ oxidoreductase', 'IDH/NADP+'), 
           c('2-oxoglutarate dehydrogenase', 'OGDH'), 
           c('Succinyl-CoA ligase (GDP-forming)', 'SCS/GTP'), 
           c('Succinate  CoA ligase  ADP forming', 'SCS/ATP'), 
           c('Succinate dehydrogenase', 'SDH'), 
           c('Fumarase  mitochondrial', 'FH mito'), 
           c('Fumarase', 'FH'), 
           c('Malate dehydrogenase  mitochondrial', 'MDH mito'), 
           c('Malate dehydrogenase', 'MDH'), 
           c('Citrate synthase', 'CS'), 
           c('ATP-Citrate lyase', 'ACLY'), 
           c('Citrate lyase', 'CLY'))
    
    for(i in replacements) {
      
      x <- stri_replace(x, fixed = i[1], replacement = i[2])
      
    }
    
    return(x)
    
  }
  
  trycats_labeller <- function(x) {
    
    replacements <- 
      list(c('3 hydroxyanthranilate 3 4 dioxygenase', 
             'HAAO'), 
           c('3-Hydroxy-L-kynurenine:2-oxoglutarate aminotransferase', 
             '3HKYNAKGAT'), 
           c('3-Hydroxykynurenamine decarboxy-lyase', 
             '3HXKYNDCL'), 
           c('3-Hydroxykynurenamine:oxygen oxidoreductase', 
             '3HXKYNOXDA'), 
           c('5-Hydroxy-L-tryptophan:oxygen 2,3-dioxygenase ', 
             'IDO'), 
           c('GTP cyclohydrolase I', 
             'GTPCH'), 
           c('3 Hydroxy L kynurenine hydrolase', 
             'KYNU'), 
           c('Kynureninase', 
             'KYNU'), 
           c('Kynurenine 3 monooxygenase', 'KMO'), 
           c('L-Kynurenine:2-oxoglutarate aminotransferase', 'KAT'), 
           c('6-pyruvoyltetrahydropterin synthase', 'PTPS'), 
           c('Sepiapterin reductase', 'SPR'), 
           c('L Tryptophanoxygen 2 3 oxidoreductase  decyclizing', 
             'TDO'), 
           c('N Formyl L kynurenine amidohydrolase', 'FKYNH'), 
           c('4--2,4-dioxobutanoate dehydratase', 'KYNATESYN'), 
           c('Quinolinate Synthase  Eukaryotic', 'QS'))
    
    for(i in replacements) {
      
      x <- stri_replace(x, fixed = i[1], replacement = i[2])
      
    }
    
    return(x)
    
  }
  
# END -----  

