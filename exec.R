# Source the complete analysis pipeline

  library(soucer)
  
  print(source_all(c('import.R', 
                     'analysis.R', 
                     'survival.R', 
                     'clustering.R', 
                     'dge.R', 
                     'metabolism.R', 
                     'paper.R'), 
                   message = TRUE, crash = TRUE))