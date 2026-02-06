library(tidyverse)
read_tsv("Orthogroups.GeneCount.tsv") %>% 
  dplyr::select(-Total) %>% 
  column_to_rownames("Orthogroup") %>% 
  mutate(across(everything(),~ifelse(.>0,1,0))) %>% 
  rowSums() %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  magrittr::set_colnames(c("familyID","total")) %>% 
  mutate(group=case_when(
    total == 17 ~ "Core",
    total < 17 & total >=15 ~ "SoftCore",
    total < 15 & total >=2 ~ "Dispensable",
    total == 1 ~ "Private"
  )) -> dat.family.group

dat.family.group %>% 
  pull(group) %>% 
  table() %>% 
  as.data.frame() %>% 
  magrittr::set_colnames(c("group","Freq")) %>% 
  mutate(group=factor(group,levels=c("Core","SoftCore",
                                     "Dispensable","Private"))) %>% 
  arrange(group) %>% 
  mutate(x=cumsum(Freq)+0.5) -> segment.df

library(ggrastr)
pdf(file = "PAV_heat_map.pdf",
    width = 12,height = 6)
read_tsv("Orthogroups.GeneCount.tsv") %>% 
  dplyr::select(-Total) %>% 
  column_to_rownames("Orthogroup") %>% 
  mutate(across(everything(),~ifelse(.>0,1,0))) %>%
  rownames_to_column("familyID") %>% 
  left_join(dat.family.group) %>% 
  arrange(desc(total)) %>% 
  mutate(x=row_number()) %>% 
  dplyr::select(-total,-group,-familyID) %>% 
  pivot_longer(!x,names_to = "sample_id",
               values_to = "group") %>% 
  mutate(group=factor(group)) %>% 
  ggplot(aes(x=x,y=sample_id))+
  geom_tile_rast(aes(fill=group))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  scale_x_continuous(expand = expansion(mult = c(0,0)))+
  scale_fill_manual(values = c("0"="#3d81c3",
                               "1"="#c25f64"),
                    labels=c("0"="Absent",
                             "1"="Present"),
                    name=NULL)+
  geom_segment(data = segment.df,
               aes(x=x,xend = x,y=-Inf,yend = Inf),
               lty="dashed",
               color="grey")+
  labs(x=NULL,y=NULL)
dev.off()
