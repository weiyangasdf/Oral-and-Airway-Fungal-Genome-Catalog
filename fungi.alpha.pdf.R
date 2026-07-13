setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/400.species.alpha")
rm(list=ls())

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}}

numFunc = function(x){
  if(x < 0.01){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else if(x<0.05){formatC(x, digits = 3, width = 1, format = "f", flag = "0")}
  else{NA}
}

numFunc_rev = function(x){
  if(x >= 0.05){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else{NA}
}


top_sites = c("Saliva","Sputum","BALF","Dental plaque","Tongue","Oropharynx", "Subgingival plaque")
load("../00.data/colors.RData")

dt = read.table("./fungi.alpha.tsv", sep=",")
samp = read.table("../00.data/sample.info", sep="\t", header=T, check.names=F )
samp = samp %>% filter(Sample_site %in% top_sites, Sample_ID %in% rownames(dt)) %>%
  group_by(BioProject_ID) %>%
  filter(n()>=10) %>%
  ungroup()


dtf = merge(samp, dt, by.x='Sample_ID', by.y='row.names')
pdt = dtf %>%
  mutate(x = paste(BioProject_ID,"(", Sample_site,")", sep=""),
         Sample_site = factor(Sample_site, levels=top_sites)) %>%
  arrange(Sample_site)

pdt$x = factor(pdt$x, levels=unique(pdt$x))

p1 <- ggplot(pdt, aes(x=x, y=shann, fill=Sample_site))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values=colors.top_sites)

p2 <- ggplot(pdt, aes(x=x, y=obs, fill=Sample_site))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values=colors.top_sites)

comps = combn(top_sites, 2, list)

p3 <- ggplot(pdt, aes(x=Sample_site, y=shann, fill=Sample_site))+
  geom_boxplot()+
  geom_signif(comparisons = comps, step_increase = 0.1, tip_length = NA,map_signif_level=numFunc_rev)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values=colors.top_sites)


p4 <- ggplot(pdt, aes(x=Sample_site, y=obs, fill=Sample_site))+
  geom_boxplot()+
  geom_signif(comparisons = comps, step_increase = 0.1, tip_length = NA,map_signif_level=numFunc_rev)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_fill_manual(values=colors.top_sites)


p <- ggpubr::ggarrange(plotlist=list(p1,p3,p2,p4), ncol=2, nrow=2,
                       widths=c(3,1))


# ggsave("./fungi.alpha.pdf", p, width=12, height=8)

