setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/201.abun.preva")
rm(list=ls())

dt = read.table("../00.data/merged.rc.norm.euk-airway.real", sep="\t", header=T, check.names=F, row.names=1)
sgbmap = read.table("../00.data/fungi.taxo-airway", sep="\t", header=T)
rownames(sgbmap) = sgbmap$species
sgbmap = sgbmap[rownames(dt),]
tmp.color = read.table("./family.color.tsv", sep="\t", header=T, comment.char = "")
colors.family = c(tmp.color$color)
names(colors.family) = tmp.color$family

pdt = data.frame(obs = rowSums(dt>0), mra = rowMeans(dt))
pdt$rate = pdt$obs / ncol(dt)
pdt$size = pdt$mra * pdt$rate
pdt$taxo = sgbmap[,'family']
pdt$label = ifelse(pdt$rate>0.5, rownames(pdt), NA)

px = ggplot(pdt, aes(x=rate))+
  geom_histogram(bins=100)+
  theme_bw()


py <- ggplot(pdt, aes(y=mra))+
  geom_histogram(bins=100)+
  scale_y_continuous(trans='log10')+
  theme_bw()


pm <- ggplot(pdt, aes(x=rate, y=mra, size=size, color=taxo))+
  geom_point(aes(size=size))+
  geom_text_repel(aes(label = label), size=4)+
  scale_y_continuous(trans='log10')+
  scale_size_continuous(trans='sqrt')+
  scale_color_manual(values=colors.family)+
  theme_bw()

p <- ggpubr::ggarrange(plotlist=list(px, ggplot(), pm, py), heights=c(1,3), widths=c(3,1), ncol=2, nrow=2)

ggsave(filename="prevalence.mra.pdf", plot=p, width=12, height=12)
