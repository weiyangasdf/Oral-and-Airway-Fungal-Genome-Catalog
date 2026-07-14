setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/202.abun.PCoA//")
rm(list=ls())

source("/share/data1/zhangy2/scripts/R_my_functions/zy_PCoA.R")
samp = read.table("../200.abun.comp/sample.to.proj", sep="\t", header=T, check.names=F)
colnames(samp) = c("sample","proj","site")

## get adonis
get_ados <- function(dt, sampf){
  
  groups = table(sampf$site)
  res = rbind()
  for(j in names(groups)){
    for(k in names(groups)){
      if(j == k){
        r2=NA
        p=NA
      }else if( (groups[[j]] < 3 || groups[[k]] < 3)){
        r2=NA
        p=NA
      }else{
        tmp.map = subset(sampf, site %in% c(j,k))
        dtf = dt[,tmp.map$proj]
        ados = adonis2(t(dtf) ~ tmp.map$group)
        r2 = ados$R2[1]
        p  = ados$`Pr(>F)`[1]
      }
      tmp.data = data.frame(a=j,b=k,pvalue=p, r2=r2)
      res = rbind(tmp.data, res)
    }
  }
  res %>% mutate(sig = ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05, "*",NA)))
}

top_sites = c("Saliva","Sputum","BALF","Dental plaque","Tongue","Oropharynx", "Subgingival plaque")
load("../00.data/colors.RData")

### genus
dt = read.table("../200.abun.comp/prof.genus.proj", sep="\t", header=T, check.names=F, row.names=1)
sampf = unique(samp[,c("proj","site")]) %>% filter(proj %in% colnames(dt))
sampf <- sampf %>% mutate(group = ifelse(site %in% top_sites, site, 'other')) %>%
  filter(group!="other")


p1 <- zy_pcoa(dt, sampf, group="group", ID="proj", title = "PCoA genus", sample.color = colors.top_sites)
p1

######### adonis
ados_g = get_ados(dt, sampf)


### species
dt = read.table("../200.abun.comp/prof.species.proj", sep="\t", header=T, check.names=F, row.names=1)
p2 <- zy_pcoa(dt, sampf, group="group", ID="proj", title="PCoA species", sample.color = colors.top_sites)
ados_s = get_ados(dt, sampf)


head(ados_s)

pdt = ados_s
p3 <- ggplot(pdt, aes(x=a,y=b,size=r2, color=r2))+
  geom_point()+
  geom_text(aes(label=sig), color="black", size=3)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_size_continuous(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        range = c(0, 10))+
  scale_color_continuous(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                         palette = c("#FEE0D2", "#FC9272", "#DE2D26", "black"))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

pdt = ados_g
p4 <- ggplot(pdt, aes(x=a,y=b,size=r2, color=r2))+
  geom_point()+
  geom_text(aes(label=sig), color="black", size=3)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_size_continuous(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                        range = c(0, 10))+
  scale_color_continuous(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                         palette = c("#FEE0D2", "#FC9272", "#DE2D26", "black"))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
p4


p <- ggpubr::ggarrange(plotlist=list(p1$plot+theme(aspect.ratio = 1),p3,
                                     p2$plot+theme(aspect.ratio = 1),p4),
                       ncol=2, nrow=2)
p

ggsave("pcoa.project.pdf", p, width=12, height=12)
