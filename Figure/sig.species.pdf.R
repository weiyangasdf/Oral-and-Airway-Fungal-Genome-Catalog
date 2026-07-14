setwd("/share/data1/zhangy2/project/79.oral_culture_fungi/aa.analyses/601.disease.meta")
rm(list=ls())

library(ComplexHeatmap)
library(circlize)

load("./maaslin.res.RData")

fm = data.frame(orogin=rownames(dtf), format = make.names(rownames(dtf)))
head(fm)
resff = merge(resf, fm, by.x='feature', by.y='format')

dtff = dtf[resff$orogin,]

sampf$name = paste(sampf$BioProject_ID, ";", sampf$Sample_site, ";", sampf$Disease_name, sep="")

res.plot = matrix(data=NA, 
                  ncol=length(unique(sampf$name)),
                  nrow = nrow(dtff),
                  dimnames=list(rownames(dtff), unique(sampf$name)))
sig.plot = res.plot

for( proj in unique(sampf$name)){
  
  x = subset( sampf, name == proj & Disease == "Control" )
  y = subset( sampf, name == proj & Disease == "Disease" )
  
  if((nrow(x)) <= 1 | (nrow(y) <= 1)){next}
  
  X = dtff[, x$Sample_ID]
  Y = dtff[, y$Sample_ID]
  
  for(j in rownames(X)){
    a = as.numeric(X[j,])
    b = as.numeric(Y[j,])
    p = wilcox.test(a,b)$p.value
    sig.plot[j,proj] = p
  }
  
  res.plot[,proj] = log2(rowMeans(Y)/rowMeans(X))
  
}

sig.plot[sig.plot<0.05] = "*"
sig.plot[sig.plot!="*"] = NA

res.plot[is.nan(res.plot)] = NA

res.plot = res.plot[, colSums(!is.na(res.plot)) > 0]

heat_matrix = res.plot
heat_matrix[heat_matrix>10] = 10
heat_matrix[heat_matrix<(-10)] = -10
mycol = colorRamp2(c(-5, -2, -1, 0, 1, 2, 5), 
           rev(c("#d73027", "#f46d43","#fee090","white","#abd9e9","#74add1","#4575b4"))
           )

nh = nrow(heat_matrix)
nw = ncol(heat_matrix)
myunit = 5
max_nchar = max(nchar(rownames(heat_matrix)))

p <- Heatmap(heat_matrix,
        
        show_row_dend = F,
        show_column_dend = F,
        na_col="white",
        border=T,
        
        rect_gp = gpar(col = "grey", lwd = 0.1), # 内部线条颜色
        height = nh * unit(myunit,"mm"), width = nw * unit(myunit,"mm"), # 保持单元格是方的
        row_names_max_width = unit(max_nchar, "char"), # 有时候legend和标签文字会重叠
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(!is.na(sig.plot[i,j])){
            grid.text(sprintf("%s", sig.plot[i, j]), x, y, gp = gpar(fontsize = 10))
          }
        },
          
        col=mycol)

pdf("sig.species.pdf", width=12, height=8)
p
dev.off()
