setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/203.abun.pvalue")
rm(list=ls())
library(Maaslin2)
library(dplyr)

dt = read.table("../00.data/merged.rc.norm.euk-airway.real.genus", sep="\t", header=T, check.names=F, row.names=1)
samp = read.table("../00.data/sample.info", sep="\t", header=T, check.names=F)

sites = c("Saliva","Sputum","BALF","Dental plaque","Tongue","Oropharynx", "Subgingival plaque")

xx = intersect(samp$Sample_ID, colnames(dt))

sampf = subset(samp, Sample_ID %in% xx & Sample_site %in% sites)
rownames(sampf) = sampf$Sample_ID

dtf = dt[, sampf$Sample_ID]
sampf = sampf %>%
  mutate(Age = ifelse(Age=="", NA,Age),
         Gender = ifelse(Gender=="",NA,Gender),
         BMI = ifelse(BMI=="",NA,BMI))

res = rbind()
for(site in sites){
  
  res.tmp <- Maaslin2(input_data = as.data.frame(t(dtf)),
         input_metadata = sampf,output = "out/all",
         plot_scatter=F,
         normalization = 'NONE',
         standardize = FALSE,
         fixed_effects  = c("Sample_site", "BMI", "Age", "Gender"),
         # reference = c("Sample_site,Dental plaque"),
         reference = c(paste("Sample_site,",site,sep="")),
         cores = 10, random_effects=c("BioProject_ID")
         )
  resf = res.tmp$results
  resf$ref = site
  res = rbind(res, resf)
}

write.table(res, "res.masslin.genus.tsv", sep=";")

