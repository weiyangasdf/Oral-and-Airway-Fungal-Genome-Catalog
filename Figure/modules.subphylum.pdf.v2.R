setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/002.genome2family.kegg")
rm(list=ls())

library(ComplexHeatmap)
library(circlize)

mod.p = read.table("./modules.pval.tsv", sep="\t")
mod.sig <- subset(mod.p, Complete < 0.05)
mod.desc <- read.table("modules.desc", sep="\t", header=T)
mod.desc = subset(mod.desc, name %in% rownames(mod.sig))

dt = read.table("../00.data/prof.func.kegg.mod_lack.airway.family", sep="\t", header=T, check.names=F, row.names=1)
forder = read.table("../00.data/family.order.csv", sep=",")
taxo = read.table("../00.data/fungi.genomes.taxo-airway", sep="\t", header=T, check.names=F, row.names=1)
taxo = data.frame(name=taxo$family, group=taxo$`phylum(subphylum)`)
taxof = taxo[match(forder$x, taxo$name),]

dtf = dt[rownames(mod.sig), taxof$name]

mycol = colorRamp2(breaks=c(2,1,0),c("white","#9ecae1", "#2171b5"))

heat_matrix <- dtf[mod.desc$name, taxof$name]
rownames(heat_matrix) = mod.desc$description

### 物种注释条带
taxo.color = read.table("./taxo.color.tsv", sep="\t", header=T, comment.char = "")
taxo.color = structure(c(taxo.color$color,"#d9d9d9"), names=c(taxo.color$taxo, "other"))

top_anno = HeatmapAnnotation(
  phylum = taxof$group,
  #phylum_name = anno_text(csp$group),
  
  annotation_height = unit(c( 0.5, 4), "cm"),
  
  col = list(
    phylum = taxo.color
  )
)


nh = nrow(heat_matrix)
nw = ncol(heat_matrix)
myunit = 3
max_nchar = max(nchar(rownames(heat_matrix)))


ha = Heatmap(heat_matrix
             ,show_row_dend = F
             ,row_split = factor(mod.desc$levelC)
             ,row_gap=unit(0,'mm')
             ,row_title_rot = 0
             ,row_names_rot = 0
             
             ,cluster_columns = F
            # ,column_split = as.factor(taxof$group)
             ,column_gap = unit(0,'mm')
             ,show_column_dend = F
             ,column_title_rot = 90
             ,cluster_column_slices = FALSE
             
             ,col = mycol
             ,border=T
             ,rect_gp = gpar(col = "grey", lwd = 0.1) # 内部线条颜色
             
             ,top_annotation = top_anno
             ,height = nh * unit(myunit,"mm")
             , width = nw * unit(myunit,"mm") # 保持单元格是方的
             ,row_names_max_width = unit(max_nchar, "char") # 有时候legend和标签文字会重叠
             
)


pdf("modules.subphylum.pdf", width=60 ,height=20)
draw(ha)
dev.off()
