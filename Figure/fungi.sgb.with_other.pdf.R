setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/103.sgb.clust.with_other")
rm(list=ls())

dt = read.table("clust.with.other.tsv", sep="\t")
head(dt)

pdt.list = list(
oafgc = unlist(subset(dt, V2!="", "V1")),
digestive = unlist(subset(dt, V3!="", "V1")),
other = unlist(subset(dt, V4!="", "V1"))
)


library(VennDiagram)

venn = venn.diagram(
  disable.logging = T,
  x = pdt.list,
  #filename = "interpro.venn.tiff",
  filename = NULL,
  fill=c("#f44336", "#8fce00", '#2ac6e3')
)

pdf("fungi.sgb.with_other.pdf", width=6, height=5)

grid.newpage()
grid.draw(venn)

dev.off()
