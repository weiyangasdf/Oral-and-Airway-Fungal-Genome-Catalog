setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/102.sgb.roc")
rm(list=ls())

dt = read.table("./roc.tsv", sep="\t")
head(dt)

pdt <- dt %>%
  group_by(V1,V2) %>%
  summarise(value=mean(V4))

p <- ggplot(pdt, aes(x=V2, y=value, color=V1))+
  geom_point()+
  geom_line()+
  theme_bw()

ggsave("sgb.num.roc.pdf", p, width=5, height=3)
