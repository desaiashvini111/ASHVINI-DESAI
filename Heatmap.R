path1 = "C:/Users/Ashvini Desai/Desktop/atg7"
setwd(path1)

m <- as.data.frame(fread("go_HEATMAP_UP_INPUT.txt"))
head(m)
row.names(m) <- m$V1
head(row.names(m))
expression = m[2:length(m[1,])]
rownames(expression) = make.names(m[,1], unique=TRUE) 
expression[1:4,1:4]
rownames=1

##install.packages(c("R.basic"), contriburl="http://www.braju.com/R/repos/")
library(R.basic)
for (i in 1:length(expression[,1])){
  expression[i,] <- zscore(as.numeric(expression[i,]))
}
length(expression[1,])
nrow(expression)

## heatmap 2
library(gplots)
mRNA_palette = greenred(64)
break.n = seq(-1.5,1.5,3/64)
set1 = as.matrix(c(rep('N', 2), rep('O', 2)))
set1[grep('N',set1)] = 5
set1[grep('O',set1)] = 2
as.factor(set1)
hc_row = hclust(dist(as.matrix(expression)),"com")
hc_col = hclust(dist(as.matrix(t(expression))),"com")
c_hc_row_adipo = cutree(hc_row, k=1)
heatmap.2(as.matrix(expression),xlab = "", ylab = "", col = mRNA_palette, trace = "none",breaks = break.n,
          Rowv=as.dendrogram(hc_row), Colv= as.dendrogram(hc_col), cexCol = 1.8, cexRow = 0.5,
          keysize = 1.0, density.info = "density", key.par=list(mgp=c(2.0,1,0), mar=c(2, 4, 4, 0.3)),
          margins = c(5,5),symm=FALSE,symkey=TRUE,symbreaks=TRUE, scale = "none", labRow = NA, colsep = NULL,
          sepwidth = 0.005)#, ColSideColors = adipo_g0)
