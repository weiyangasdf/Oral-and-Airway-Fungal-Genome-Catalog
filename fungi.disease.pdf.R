setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/600.disease")
rm(list=ls())

load("./plot.RData")

f1 = samp %>% dplyr::select(Disease_name, BioProject_ID) %>% unique() %>% group_by(Disease_name) %>% summarise(n = n()) %>% arrange(desc(n))

samp = samp %>%
  mutate(Disease_name = factor(Disease_name, levels=f1$Disease_name)) %>%
  arrange(desc(Disease_name), desc(Sample_site))


colors.disease = c("#45bcec", "#f591b6","#efd493", "#e44a33", "#8c83b4", "#bfef45","#000075","#cdaa7d", "#ffcdcc", "#800000", "#bebedb", "#34a853")
names(colors.disease) = c("Dental_caries", "Periodontitis", "Upper_respiratory_tract_infection", "COVID","Pneumonia",
                          "Chronic obstructive pulmonary disease", "Cystic_fibrosis", "Pancreatic_ductal_carcinoma",
                          "Colorectal_cancer", "Rheumatoid_arthritis", "Type_1_diabetes_mellitu", "Autism_spectrum_disorder")

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}


## obs

pdt = merge(res.alpha, unique(samp[,c("name","Disease_name")]), by.x='name', all.x=T)
pdt$name = factor(pdt$name, levels=unique(samp$name))
pdt$sig = unlist(lapply(pdt$obs.p, sigFunc))

p1 <- ggplot(pdt, aes(x=obs, y=name, color=Disease_name))+
  geom_segment(aes(y=name, x=0,yend=name, xend=obs))+
  geom_point(size=2)+
  geom_text(aes(label=sig), color="red")+
  scale_color_manual(values=colors.disease)+
  labs(x="OBS", y=NULL)


## shann
pdt$sig = unlist(lapply(pdt$shann.p, sigFunc))
p2 <- ggplot(pdt, aes(x=shann, y=name, color=Disease_name))+
  geom_segment(aes(y=name, x=0,yend=name, xend=shann))+
  geom_point(size=2)+
  geom_text(aes(label=sig), color="red")+
  scale_color_manual(values=colors.disease)+
  labs(x="Shann", y=NULL)

## R2
pdt = merge(res.r2, unique(samp[,c("name","Disease_name")]), by.x='name', all.x=T)
pdt$name = factor(pdt$name, levels=unique(samp$name))
pdt$sig = unlist(lapply(pdt$p, sigFunc))

p3 <- ggplot(pdt, aes(x=r2, y=name, color=Disease_name))+
  geom_segment(aes(y=name, x=0,yend=name, xend=r2))+
  geom_point(size=2)+
  geom_text(aes(label=sig), color="red")+
  scale_color_manual(values=colors.disease)+
  labs(x="R2", y=NULL)

## auc
pdt = merge(aucs$table, unique(samp[,c("name","Disease_name")]), by.x='name', all.x=T)
pdt$name = factor(pdt$name, levels=unique(samp$name))
p4 <- ggplot(pdt, aes(x=auc, y=name, fill=Disease_name))+
  geom_bar(stat='identity')+
  geom_errorbar(aes(xmin=low, xmax=high), width=NA)+
  geom_vline(xintercept = c(50, 70), lty="dashed")+
  scale_fill_manual(values=colors.disease)+
  labs(x="AUC", y=NULL)


ot = theme(axis.text.y=element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = 'none')

p <- ggpubr::ggarrange(plotlist=list(p1,p1+ot,p2+ot,p3+ot,p4+ot), ncol=5)
p

ggsave("fungi.disease.pdf", p, width=15, height=8)