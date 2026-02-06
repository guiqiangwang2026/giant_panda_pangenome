library(tidyverse)
nodemat<-"./nodemat.tsv"
datmat <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)

datmat$combres <- rowSums(datmat %>% select(-nodeid))
datmat
totassemb <- ncol(datmat) - 2
totassemb

graphlen<-"./graph_len.tsv"
datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")
datlen
datmat <- datmat %>% left_join(datlen, by = c("nodeid"))
datmat %>% head()

core <- sum(datmat[datmat$combres == totassemb, "conlen"])
core


datexp <- read.table(nodemat, header = TRUE, stringsAsFactors = FALSE)
datexp
datlen <- read.table(graphlen, header = FALSE, stringsAsFactors = FALSE)
datlen %>% head()

colnames(datlen) <- c("nodeid", "conlen", "chromo", "pos", "rrank")


breeds <- colnames(datexp %>% select(-nodeid))
breeds
no_rep <- ncol(datexp) - 1
no_rep


datpan <- data.frame(
  norep = numeric(),
  nosamp = numeric(),
  assemb = character(),
  core_gen = numeric(),
  flex_gen = numeric(),
  tot_gen = numeric()
)

datcon <- datlen %>% select(nodeid, conlen)
datcon


for (nosamp in seq(1, length(breeds))) {
  # how many sampling repeated
  for (norep in seq(1, no_rep)) {
    selsamp <- as.character(sample(breeds, size = nosamp))
    # selected breeds
    selcol <- datexp[, colnames(datexp) %in% selsamp, drop = FALSE]
    # add contig length
    # it is ordered so we just add it from conlen
    # add number of colour in node
    selcol$comcol <- rowSums(selcol)
    # add contig len
    selcol$conlen <- datlen$conlen
    # core genome
    # core if shared in all of the member of population
    core_gen <- selcol[selcol$comcol == nosamp, "conlen"] %>% sum()
    # flex genome if shared less than the total of population
    # not consider node not present 
    if (nosamp == 1) 
      flex_gen  <- 0 
    else 
      flex_gen <- selcol[selcol$comcol < nosamp & selcol$comcol > 0, "conlen"] %>% sum()
    # tot_genome if at least observed in a single breeds
    tot_gen <- selcol[selcol$comcol > 0, "conlen"] %>% sum()
    datpan <- rbind(datpan, data.frame(
      norep = norep,
      nosamp = nosamp,
      assemb = paste(selsamp, collapse = ","),
      core_gen = core_gen,
      flex_gen = flex_gen,
      tot_gen = tot_gen
    ))
  }
}

datpan


p1<-datpan %>% 
  select(nosamp,core_gen,tot_gen) %>% 
  mutate(nosamp=factor(nosamp)) %>% 
  pivot_longer(!nosamp) %>% 
  ggplot(aes(x=nosamp,y=value))+
  geom_boxplot(aes(fill=name))+
  geom_smooth(aes(x=as.numeric(nosamp),color=name))+
  #geom_violin(aes(fill=name))+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank())+
  scale_y_continuous(labels = function(x){paste(x/1000000,"M")})+
  labs(x="Sample Number",y="Genome Size")+
  scale_color_manual(values = c("#e3010a","#00b2ec"))+
  scale_fill_manual(values = c("#e3010a","#00b2ec"))

p2<-datpan %>% 
  select(nosamp,core_gen,tot_gen) %>% 
  mutate(nosamp=factor(nosamp)) %>% 
  pivot_longer(!nosamp) %>% 
  ggplot(aes(x=nosamp,y=value,color=name))+
  geom_jitter(width = 0.1,size=5,alpha=0.8)+
  geom_smooth(aes(x=as.numeric(nosamp)))+
  #geom_violin(aes(fill=name))+
  theme_bw(base_size = 20)+
  theme(panel.grid = element_blank(),
        legend.position = "top",
        legend.title = element_blank())+
  scale_y_continuous(labels = function(x){paste(x/1000000,"M")})+
  labs(x="Sample Number",y="Genome Size")+
  scale_fill_manual(values = c("#4da0a0","#9b3a74"))+
  scale_color_manual(values=c("#4da0a0","#9b3a74"))

library(patchwork)
p1+p2
