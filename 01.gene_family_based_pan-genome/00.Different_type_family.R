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
# 导出基因家族分类文件到本地
# install.packages("writexl")
library(writexl)
write_xlsx(dat.family.group, "different.family.xlsx")

myfun<-function(x){
  y<-c()
  for(i in 1:length(x)){
    if(i == 1){
      new_y = x[i]/2
      y<-append(y,new_y)
    }
    else{
      new_y = sum(x[1:i-1])+x[i]/2
      y<-append(y,new_y)
    }
  }
  return(y)
}

myfun(c(1,4,7,8))

library(ggrepel)

dat.family.group %>% 
  group_by(group) %>% 
  summarise(count=n()) %>% 
  mutate(prop=sprintf("%0.2f",count*100/sum(count))) %>% 
  mutate(label=paste0(prop,"%")) %>% 
  mutate(group=factor(group,levels=rev(c("Core","SoftCore","Dispensable","Private")))) %>% 
  arrange(group) %>% 
  ungroup() %>% 
  mutate(y.position=myfun(rev(count))) %>% 
  ggplot(aes(y=count,x="",fill=group))+
  geom_text_repel(aes(label=paste0(rev(group),
                                   "\n",
                                   rev(label)),
                      y=y.position),
                  nudge_x = c(1,1,2,1))+
  geom_bar(width = 1, stat="identity", alpha=1,
           show.legend = FALSE)+
  coord_polar("y", start=0)+
  theme_void(base_size = 15)+
  scale_fill_manual(values = rev(c("#e46a6f","#e58ea2","#5f9ed5","#e6db91")))+
  theme(plot.margin = unit(c(1,1,1,1),'cm'))-> fig2.1

fig2.1

data <- read.table("Orthogroups_UnassignedGenes.tsv", header = TRUE, sep = "\t")
no_first_column <- data[-1]
non_empty <- apply(no_first_column, 2, function(x) sum(!is.na(x) & x != ""))
bar.dat.02 <- data.frame(species = names(no_first_column), Freq = non_empty)

dat.family.group %>% 
  pull(total) %>% table() %>% 
  as.data.frame() %>% 
  magrittr::set_colnames(c("Frequency","Family number")) %>% 
  mutate(Frequency=as.numeric(Frequency)) %>% 
  mutate(group=case_when(
    Frequency == 17 ~ "Core",
    Frequency < 17 & Frequency >=15 ~ "SoftCore",
    Frequency < 15 & Frequency >=2 ~ "Dispensable",
    Frequency == 1 ~ "Private"
  )) %>% 
  mutate(group=factor(group,levels=c("Core","SoftCore","Dispensable","Private"))) -> bar.dat.01

ggplot(data=bar.dat.01,
       aes(x=Frequency,y=`Family number`))+
  geom_col(aes(fill = group),
           show.legend = FALSE)+
  scale_x_continuous(breaks = c(seq(1,27,3),27))+
  scale_fill_manual(values = c("#e46a6f","#e58ea2","#5f9ed5","#e6db91"))+
  theme_bw(base_size = 15)+
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), breaks = seq(0, 20000, 500))+
  geom_bar(data=bar.dat.02,
           aes(x=1,y=Freq),
           stat="identity",
           position = "stack",
           color="grey",
           lty="dashed",
           show.legend = FALSE,
           fill="#e6db91")-> fig2.2

fig2.2
