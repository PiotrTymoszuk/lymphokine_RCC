# Recon2 metabolic and Reactome signatures: calculation of ssGSEA scores

  insert_head()
  
# parallel backend ------
  
  insert_msg('Parallel backend')
  
  plan('multisession')
  
# expression data frames, log scale ------
  
  insert_msg('Preparing expression data frames')
  
  sig$annotation <- list(tcga = tcga$annotation$gene_symbol,
                         emtab1980 = emtab1980$annotation$gene_symbol,
                         gse73731 = gse73731$annotation$gene_symbol,
                         gse167093 = gse167093$annotation$gene_symbol,
                         reca = reca$annotation$gene_symbol,
                         cm10 = cm10$annotation$gene_symbol,
                         cm25ev = cm25ev$annotation$gene_symbol,
                         cm25ni = cm25ni$annotation$gene_symbol)
  
  sig$expression <- list(tcga = tcga$expression %>%
                           filter(tissue_type == 'Tumor'),
                         emtab1980 = emtab1980$expression,
                         gse73731 = gse73731$expression,
                         gse167093 = gse167093$expression %>%
                           filter(tissue_type == 'Tumor'),
                         reca = reca$expression,
                         cm10 = cm10$expression,
                         cm25ev = cm25ev$expression,
                         cm25ni = cm25ni$expression)
  
  sig$expression <- map2(sig$expression,
                            sig$annotation,
                            ~select(.x, patient_id, any_of(.y)))
  
  sig$patient_id <- sig$expression %>% 
    map(~.x$patient_id)
  
  sig$expression <- sig$expression %>% 
    map(select, - patient_id)
  
# Reactome signatures -------
  
  insert_msg('Reactome signatures')
  
  ## loading the local MSig database
  
  sig$reactome_db <- 
    load_dbsig('./data/signatures/msigdb.v7.5.1.symbols.gmt') %>% 
    search(key = 'sign_name', 
           value = '^REACTOME_')
  
  ## calculation, calling the method directly

  sig$reactome_scores <- sig$expression %>% 
    future_map(~calculate.dbsig(x = sig$reactome_db, data = .x))
  
  sig$reactome_scores <- 
    map2(sig$reactome_scores, 
         sig$patient_id, 
         ~mutate(.x, 
                 patient_id = .y, 
                 .before = .data[[names(.x)[1]]]))
  
  ## signature lexicon
  
  sig$reactome_lexicon <- sig$reactome_db$sign_name %>% 
    stri_replace(regex = '^REACTOME_', replacement = '') %>% 
    stri_replace_all(fixed = '_', replacement = ' ') %>% 
    set_names(sig$reactome_db$sign_name) %>% 
    compress(names_to = 'variable', 
             values_to = 'label')
  
# Reacon signatures ------
  
  insert_msg('Recon signatures')
  
  ## extracting genes associated with the Recon subsystems
  
  sig$recon_lexicon$subs <- extract_subsystems(Recon2D)
  
  sig$recon_lexicon$genes <- extract_genes(Recon2D)
  
  ## not all reactions have genes assigned to
  ## hence genes on the LHS
  
  sig$recon_lexicon <- left_join(sig$recon_lexicon$genes, 
                                 sig$recon_lexicon$subs, by = 'react_id')
  
  sig$recon_lexicon <- sig$recon_lexicon %>% 
    mutate(variable = make.names(subsystem), 
           variable = tolower(variable), 
           variable = stri_replace_all(variable, 
                                       regex = '\\.+', 
                                       replacement = '_'), 
           label = subsystem, 
           gene_symbol = exchange(entrez_id, 
                                  dict = tcga$annotation, 
                                  key = 'entrez_id', 
                                  value = 'gene_symbol'))
  
  ## gene list
  
  sig$recon_db <- sig$recon_lexicon %>% 
    dlply('variable', function(x) x$gene_symbol) %>% 
    map(unname) %>% 
    map(reduce, union)
  
  ## calculation of the signatures
  
  sig$recon_scores <- sig$expression %>% 
    future_map(~calculate.default(x = sig$recon_db, data = .x))
  
  sig$recon_scores <- 
    map2(sig$recon_scores, 
         sig$patient_id, 
         ~mutate(.x, 
                 patient_id = .y, 
                 .before = .data[[names(.x)[1]]]))
  
# catching the signture scores and their lexicons -------
  
  insert_msg('Caching the signature data')
  
  sig <- sig[c('reactome_scores', 
               'reactome_lexicon', 
               'recon_scores', 
               'recon_lexicon')]
  
  save(sig, file = './data/signatures/signatures.RData')
    
# END ------
  
  plan('sequential')
  
  insert_tail()