setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/001.genomes.info")
rm(list=ls())

dt = read.table("../00.data/fungi.genomes.info", sep="\t", header=T)
dt2 = dt
dt2$source = "all"
pdt = rbind(dt, dt2)

colors.source = c("#f5f5b3","#bab6d5","#89cabe","gray")
levs = c("Isolated(pub)","Isolated(this study)", "MAG(this study)", "all")
names(colors.source) = levs

pdt$source = factor(pdt$source, levels=levs)

## Genomesize
p1 <- ggplot(pdt, aes(x=source, y=sum_len, fill=source ))+
  geom_violin()+
  geom_boxplot(width=0.05, outlier.color=NA, fill="white")+
  scale_y_continuous(n.breaks = 10)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=NULL, y="Genome Size(bp)")+
  scale_fill_manual(values=colors.source)
p1

## N50
p2 <- ggplot(pdt, aes(x=source, y=N50, fill=source ))+
  geom_violin()+
  geom_boxplot(width=0.05, outlier.color=NA, fill="white")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=NULL, y="N50 length(bp)")+
  scale_fill_manual(values=colors.source)



## GC
p3 <- ggplot(pdt, aes(x=source, y=GC, fill=source ))+
  geom_violin()+
  geom_boxplot(width=0.05, outlier.color=NA, fill="white")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x=NULL, y="(%)GC")+
  scale_fill_manual(values=colors.source)


p <- ggpubr::ggarrange(plotlist=list(p1,p2,p3), ncol=3, nrow=1)
ggsave("./fungi.info.pdf", p, width=15, height=3)
