setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/004.genome.eggnog")
rm(list=ls())

dt = read.table("../00.data/prof.func.eggnog.category", sep="\t", header=T, check.names=F, row.names=1)
sgb.map = read.table("../00.data/fungi.genomes.taxo-airway", sep="\t", header=T, check.names=F)
egg.class = read.table("/share/data1/Database/eggnog/eggnog.func.taxo", sep="\t")
taxo.color = read.table("../003.genome2species.kegg/taxo.color.tsv", sep="\t", comment.char = "", header=T)
color.taxo = structure(taxo.color$color, names=taxo.color$taxo)

dtf = t(dt[,sgb.map$file])

taxof = table(sgb.map$`phylum(subphylum)`)
taxof = names(taxof[taxof>=10])

pdt <- aggregate(dtf, by=list(sgb.map$`phylum(subphylum)`), mean) %>%
  melt() %>%
  merge(., egg.class, by.x='variable', by.y='V1', all.x=T) %>%
  filter(!(variable %in% c("-","S")), Group.1 %in% taxof) %>%
  mutate(V2 = factor(V2, levels=c("INFORMATION STORAGE AND PROCESSING","CELLULAR PROCESSES AND SIGNALING","METABOLISM"))) %>%
  arrange(as.numeric(V2), variable)
pdt$variable = factor(pdt$variable, levels=unique(pdt$variable))

p <- ggplot(pdt, aes(x=variable, y=value, fill=Group.1))+
  geom_bar(stat="identity")+
  facet_grid(Group.1 ~ ., scales="free_y")+
  theme_bw()+
  scale_fill_manual(values=color.taxo)+
  theme(panel.grid = element_blank())+
  labs(x=NULL, y="Numboer of orthologous groups")

ggsave("fungi.egg_catagory.pdf", p, width=7, height=4)
