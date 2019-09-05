if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("pd.hg.u95av2", version = "3.8")
library("'pd.hg.u95av2")

BiocManager::install("oligo", version = "3.8")
library("oligo")

BiocManager::install("RankProd", version = "3.8")
library("RankProd")

BiocManager::install("hugene10sttranscriptcluster.db", version = "3.8")
library("hugene10sttranscriptcluster.db")
library("proxy")
library("methods")
library("S4Vectors")
library("Hmisc")
setwd("J:\\MICROARRAY_Analysis_Data")

##Read .CEL files
Files <- list.celfiles()
rawData <- read.celfiles(Files)
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".CEL", "", filename)
sampleNames(rawData) <- sampleNames

chip <- gsub('\\.','',rawData@annotation)
chip <- substr(chip,3,nchar(chip))
Database <- paste(chip,"transcriptcluster.db",sep="")

####Gene PreProcessing 'RMA'
GeneData <- rma(rawData)
my_data <- exprs(GeneData)

#Format values to 5 decimal places
rma <- format(my_data, digits=5)

#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
Symbols = unlist(mget(probes, hugene10sttranscriptclusterSYMBOL, 
                      ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hugene10sttranscriptclusterENTREZID, 
                         ifnotfound=NA))

#Combine gene annotations with raw data
rma = cbind(probes,Symbols,Entrez_IDs,rma)

#RMA-normalized, mapped data to file
write.table(rma, file = "RMA_Results.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE, col.names = TRUE)

require(ggplot2)
## Boxplot before preprocessing
tata <- as.data.frame(exprs(rawData))
tiff(file="Boxplot_before.tiff", bg="transparent", width=500, height=500)
boxplot2(tata, main="Box_Plot_before_Normalization",las=2, ylim=c(0,500), 
         sizeFactors=10)
dev.off()

## Boxplot after processing
tiff(file="Boxplot_after.tiff", bg="transparent", width=500, height=500)
boxplot(GeneData, main="Box_Plot_after_Normalization",las=2, sizeFactors=10)
dev.off()

Expression_set <- my_data
write.csv2(Expression_set,"ExpressionSet.csv")

# design matrix
cl <- c(rep(0,2),rep(1,2))
# ttest & fold change
source("ttest_FC.R")
FTT<- ttest_FC(Expression_set)

#Highlight genes that have an absolute fold change > 2 and a p-value < 0.05
Bonferroni cut-off
FTT$threshold = abs(FTT$FC)>0.5 & FTT$p.value<0.05


##Construct the Volcano plot object
tiff(file="Volcano_Plot.tiff", bg="transparent", width=800, height=600)
ggplot(data=FTT, aes(x=FC, y=-log10(p.value), colour=threshold)) +
  geom_point(alpha=0.6, size=2) +
  theme(legend.position = "none", axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold")) +
  xlim(c(-2, 2)) + ylim(c(0, 5)) +
  xlab("log2 fold change") + ylab("-log10 p-value")
dev.off()

#OR

##Construct the Volcano plot object
install.packages("gdata")
install.packages("gplots")
library(gdata)
library(gplots)
png(filename = "Volcanoplot.png")
with(FTT, plot(FC, -log10(p.value), pch=20, main="Volcano plot", xlim=c(-1.5,1.5), ylim=c(0, 5)))
with(subset(FTT, p.value <.05 & FC > 0.5), points(FC, -log10(p.value), pch=20, col="red"))
with(subset(FTT, p.value <.05 & FC < -0.5), points(FC, -log10(p.value), pch=20, col="green"))
dev.off()