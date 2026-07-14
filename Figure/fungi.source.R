setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/001.genomes.info")
rm(list=ls())

library(ggplot2)
library(dplyr)

source("/share/data1/zhangy2/scripts/R_my_functions/zy_pie.R")
dt = read.table("../00.data/fungi.genomes.info", sep="\t", header=T, check.names=F)
head(dt)



dtf = dt[grep("Respiratory", dt$site),]

pdt1 <- dtf %>% group_by(source) %>% count()

p1 <- zy_pie(pdt1, value="n", fill="source")
p1

p2 <- ggplot(dtf, aes(x=sum_len/1000000, y=Complete, color=source_1))+
  geom_point()+
  theme_bw()+
  scale_x_continuous()+
  labs(x="Genome Size(Mbp)", y="(%) Complete")+
  theme(panel.grid = element_blank(),
        aspect.ratio = 1)

ggplot(dtf, aes(x=sum_len/1000000,fill=source_1))+
  geom_histogram(bins=100)+
  facet_grid(source_1 ~ .)+
  theme_bw()+
  scale_x_continuous(limits=c(0,80), breaks=seq(0,100,10))+
  labs(x="Genome Size(Mbp)", y="Number of Genomes")+
  theme()


# ggsave("fungi.source.pdf", p1, width=6, height=6)
