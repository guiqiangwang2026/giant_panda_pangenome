library(tidyverse)
# read_tsv("Orthogroups.GeneCount.tsv") %>% 
  # dplyr::select(-Total) %>% 
  # column_to_rownames("Orthogroup") %>% 
  # mutate(across(everything(),~ifelse(.>0,1,0))) %>% 
  # write_delim(file = "PanGP.input",
              # delim = "",
              # col_names = FALSE)

library(ggtrendline)

PanGP.dat<-read_tsv("PanGenomeData.txt")
PanGP.dat %>% colnames()

sample.size.x<-PanGP.dat %>% pull(GenomeNum)

pan.y<-PanGP.dat %>% pull(`Pan-Genome Size`)
core.y<-PanGP.dat %>% pull(`Core Genome Size`)

ggtrendline(sample.size.x,pan.y,
            model = "exp3P",eSize = 5,
            eq.x = 20,
            eq.y = 52000,
            rrp.x = 20,
            rrp.y = 50000)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Sample size",y="Pan gene family") -> p1

ggtrendline(sample.size.x,core.y,
            model = "exp3P",eSize = 5,
            eq.x = 20,
            eq.y = 32000,
            rrp.x = 20,
            rrp.y = 30000)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Sample size",y="Core gene family") -> p2

library(patchwork)

p1+p2