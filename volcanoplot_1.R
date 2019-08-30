path1="C:/Users/Ashvini Desai/Desktop/atg7/Quantile Normalization and Volcano plot/Input file"
setwd(path1)

#install.packages("gdata")
#install.packages("gplots")
library(gdata)
library(gplots)
library(readxl)
data <- read_excel("volcano_input.xlsx", sheet=1)
#View(data)

head (data)
df <-as.data.frame(data)

png(filename = "VolcanoPlot.png")
with(df, plot(log2_FoldChange, -log10(p_value), pch=20, main="Volcano plot", xlim=c(-4,4)))
with(subset(df, p_value < .01 & log2_FoldChange > 2), points(log2_FoldChange, -log10(p_value), pch=20, col="red"))
with(subset(df, p_value < .01 & log2_FoldChange < -2), points(log2_FoldChange, -log10(p_value), pch=20, col="green"))
dev.off()
