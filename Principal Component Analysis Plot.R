m <- as.data.frame(fread("INPUT.txt"))
expression = m[2:length(m[1,])]
rownames(expression) = make.names(m[,1], unique=TRUE) 

#PCA polt ?????????
library(FactoMineR)
pca <- expression
result_pca = PCA(t(pca), graph = FALSE)
temp_adipo = pca$ind$coord
x_min_adipo = min(temp_adipo[,1]) ; x_max_adipo = max(temp_adipo[,1]) ; y_min_adipo = min(temp_adipo[,2]) ; y_max_adipo = max(temp_adipo[,2]) 
HC_adipo = temp_adipo[1:length(pca[1,]),1:2]
plot(pca)##, pch = 20)
points(HC_adipo, xlim = c(x_min_adipo-1,x_max_adipo+1), ylim = c(y_min_adipo-1,y_max_adipo+1), 
       col = c(rep(5,8), rep(2,8)), pch = 20, cex = 4)
#points(HC1, xlim = c(x_min1-1,x_max1+1), ylim = c(y_min1-1,y_max1+1), col= g0,pch = 20, cex = 4)
text(HC_adipo, rownames(HC_adipo), col = 1, cex = 1.5)

##
expression5_1 <- pairs(data.matrix(expression5_1))
