setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/006.genomes.func")
rm(list=ls())

dt = read.table("./func.anno.tsv", sep="\t", header=T)

xx = dt[,-1]/dt[,1]*100
pdt = melt(xx[,-3])

pdt$variable = factor(pdt$variable, levels=c("eggnog", "cazy","proteases","lipase","pfam","kegg"))
p <- ggplot(pdt, aes(x=variable, y=value, fill=variable))+
  geom_boxplot()+
  theme_bw()

ggsave(filename="func.annoted.pdf", plot=p, width=6, height=3)
