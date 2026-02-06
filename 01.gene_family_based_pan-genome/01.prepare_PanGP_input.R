library(tidyverse)
read_tsv("Orthogroups.GeneCount.tsv") %>% 
  dplyr::select(-Total) %>% 
  column_to_rownames("Orthogroup") %>% 
  mutate(across(everything(),~ifelse(.>0,1,0))) %>% 
  write_delim(file = "PanGP.input",
              delim = "",
              col_names = FALSE)
