setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/101.sgb.tree")
rm(list=ls())

library(picante)
tree <- read.tree("./00.tree.nwk")
map = read.table("../00.data/fungi.taxo-airway", sep="\t", header=T,check.names=F)

colors.source = c("#f5f5b3","#bab6d5","#89cabe","gray")
levs = c("new species","known(other)", "known(oral)", "all")
names(colors.source) = levs

head(map)
tmp.func <- function(taxo){

  x_count = data.frame(table(map[,taxo]))
  x_count$Var1 = as.character(x_count$Var1)
  map$xx = map[,taxo]
  map <- merge(x_count, map, by.x='Var1', by.y='xx', all.y=T) %>%
    mutate(xx = ifelse(Freq < 5, "u__other", Var1),
           xx = ifelse(xx %in% c("p__","c__","o__","f__","g__"), "u__other", xx))
  
  
  lab.ord1 = map %>%
    group_by(xx) %>%
    count()
  
  lab.ord2 = map %>%
    group_by(xx, group1) %>%
    count() %>%
    merge(lab.ord1, by="xx") %>%
    mutate(rate = n.x/n.y) %>%
    dcast(xx ~ group1, value.var="rate") %>%
    arrange(desc(`new species`), desc(`known(other)`),desc(`known(oral)`))
  
  lab.ord <- map %>%
    group_by(xx) %>%
    count() %>%
    ungroup() %>%
    mutate(label = paste(xx, "(",n,")", sep=""),
           label = gsub("^.__", "", label, perl=T),
           xx = factor(xx, levels=lab.ord2$xx)) %>%
    arrange(xx)
  
  pdt <- map %>%
    group_by(xx, group1) %>%
    count()
  
  res.culture <- merge(pdt, lab.ord, by="xx") %>%
    mutate(label = factor(label, levels=rev(lab.ord$label)),
           rate = n.x/n.y * 100)

  ################
  #     PD
  res = rbind()
  grps = unique(map$xx)
  for(grp in grps){
    x <- subset(map, xx == grp)
    y = data.frame(all=1, obj = (x$group1  == "new species")+0, row.names=x$Represent_name)
    
    pds = pd(t(y), tree, include.root=F)
    tmp <- data.frame(all_pd=pds['all',1], obj_pd = pds['obj',1], name=grp, obj_n = sum(y$obj), all_n = sum(y$all))
    res <- rbind(res, tmp)  
  }
  
  res.pd <- res %>%
    mutate(rate_pd = obj_pd / all_pd * 100,
           label = paste(name, "(",all_n,")", sep=""),
           label = gsub("^.__", "", label, perl=T)
    )
  
  pdt <- rbind(data.frame(value = res.pd$rate_pd, label=res.pd$label, all_n = res.pd$all_n, obj_n = res.pd$obj_n, fill="new species", facet="PD"),
               data.frame(value = res.culture$rate, label=res.culture$label, all_n = res.culture$n.y, obj_n = res.culture$n.x, fill=res.culture$group1, facet="% species"))
  
  pdt$label = factor(pdt$label, levels = rev(lab.ord$label))
  pdt <- pdt %>%
    mutate(text = ifelse(fill=="new species", paste(signif(value, digits = 4),"% ,", obj_n, sep=""), NA)) %>%
    mutate(rm=ifelse(fill=="new species" & facet=="PD" & obj_n<3,"Y","N")) %>%
    filter(rm=="N")
    # filter(label=="Xanthomonadales(37)")
  pdt$taxo = taxo
  
  p <- ggplot(pdt %>% filter(all_n >= 5), aes(x=value, y=label, fill=fill))+
    geom_bar(stat="identity", color="black", width=1)+
    geom_text(aes(label=text, x = 0), hjust=0)+
    facet_grid(.~facet)+
    theme_bw()+
    labs(x=NULL, y=NULL, title=taxo)+
    scale_fill_manual(values=colors.source)
  
  list(plot=p, data=pdt)
}



xx <- list(tmp.func("phylum"),
              tmp.func("class"),
              tmp.func("order"),
              tmp.func("family"),
              tmp.func("genus"))


res = rbind()

for(x in xx){
  res = rbind(res,x$data)
}

p <- ggplot(res, aes(x=value ,y = label,  fill=fill))+
  geom_bar(stat="identity", color="black", width=0.6)+
  facet_grid(taxo ~ facet, scales="free_y", space = "free_y")+
  scale_fill_manual(values=colors.source)+
  theme(axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank())
p

# ggsave("diff_taxo.PD.pdf", p, width=8, height=10)


write.table(res, "tree.pdf.tsv", sep="\t")

res %>% filter(facet=="PD") %>%
  filter(!grepl("other", label)) %>%
  group_by(taxo) %>%
  arrange(desc(value)) %>%
  filter(row_number() <= 5) %>%
  ungroup() %>%
  arrange(taxo) %>%
  print(n=50)


res %>% filter(facet=="PD") %>%
  filter(!grepl("other", label)) %>%
  group_by(taxo) %>%
  summarise(mean = mean(value), min = min(value), max = max(value) )
