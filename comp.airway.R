setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/200.abun.comp/")
rm(list=ls())


source("/share/data1/zhangy2/scripts/R_my_functions/zy_compositions.R")
samp = read.table("../00.data/sample.info", sep="\t", header=T, check.names=F)

### subphylum
dt = read.table("../00.data/merged.rc.norm.euk-airway.real.subphylum", sep="\t", header=T, check.names=F, row.names=1)
p1 <- zy_group_compositions(dt, sample_map = samp, ID="Sample_ID", group="Sample_site", top_N = 15, order_func = "order")
p1 <- p1 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("comp.airway-subphylum.pdf", p1, width=15, height=6)


### genus
dt = read.table("../00.data/merged.rc.norm.euk-airway.real.genus", sep="\t", header=T, check.names=F, row.names=1)
p2 <- zy_group_compositions(dt, sample_map = samp, ID="Sample_ID", group="Sample_site", top_N = 15, order_func = "cluster")
p2 <- p2 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("comp.airway-genus.pdf", p2, width=15, height=6)



################################################
#                   BWPR
################################################

### genus
dt = read.table("../00.data/merged.rc.norm.euk-airway.real.genus", sep="\t", header=T, check.names=F, row.names=1)
p2 <- zy_group_compositions(dt, sample_map = samp[grep("in-house|this" ,samp$NCBI_BioProject_ID),] , ID="Sample_ID", group="Sample_site", top_N = 15, order_func = "cluster")
p2 <- p2 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("comp.airway-genus-BWPR.pdf", p2, width=15, height=6)

