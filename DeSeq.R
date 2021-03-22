library("DESeq2")

mycurrentdirectory = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(mycurrentdirectory)

sampleInfo <- read.csv("shortRNAseq.txt")

sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)

## I am assuming feature counts finished
countdata <- read.delim("human_counts_RNA.txt", comment.char="#",row.names = 1)
# Remove first five columns (chr, start, end, strand, length)

countdata = countdata[ ,6:ncol(countdata)]
# Remove crap from colnames in countdata
temp = colnames(countdata)
temp = gsub("Bam.","",temp)
temp = gsub(".bam","",temp)
colnames(countdata) = temp

##  does everything match up...
cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)

# create DEseq2 object & run DEseq
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )
res <- results(dds,
               contrast = c('TissueCode','LK','LD')) #specifies which dataset to compare
res
#Had to remove seurat and sctransform: 
#remove.packages("sctransform", lib="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
#remove.packages("Seurat", lib="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")


plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )

###  throw out lowly expressed genes?? ... I leave this as an exercise
###  add external annotation to "gene ids"
# log transform
rld = rlog( dds )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)

sampleDists = dist( t( assay(rld) ) )

#install R 3.5 to avoid maybe the seurat and deseq2 conflict
# heat map
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)

#turn off developer graphics
dev.off()

# PCs
# wow you can sure tell tissue apart
print( plotPCA( rld, intgroup = "TissueCode") )
# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 400 )

png(filename='LamHeatMap.png', width=800, height=750)
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
# volcano plot this is an exercise
graphics.off()

png(filename='LKLDVolcano.png', width=800, height=750)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')
graphics.off()

###########################Volcano Plots

library(EnhancedVolcano)
library(org.Hs.eg.db)
library(dplyr)

volcanodata<-countdata
volcanodata<-volcanodata %>% mutate(Total = rowSums(.))
ens <- rownames(volcanodata)
symbols <- mapIds(org.Hs.eg.db, keys = ens,
                  column = c('SYMBOL'), keytype = 'ENSEMBL')
symbols <- symbols[!is.na(symbols)]
symbols <- symbols[match(rownames(volcanodata), names(symbols))]
combined<-cbind(symbols,volcanodata)
combined2<-combined[!is.na(combined$symbols),]
combined3<-combined2[order(combined2$Total,decreasing = TRUE),]
combined4<-combined3[!duplicated(combined3$symbols),]
rownames(combined4) <- combined4$symbols
combined4$symbols <- NULL
combined4$Total <- NULL

# create DEseq2 object & run DEseq
ddsv = DESeqDataSetFromMatrix(countData=combined4, colData=sampleInfo, design=~TissueCode)
ddsv <- DESeq(ddsv)
resv <- results( ddsv )
resv <- results(ddsv,
               contrast = c('TissueCode','LK','LD')) #specifies which dataset to compare
resv
#Had to remove seurat and sctransform: 
#remove.packages("sctransform", lib="/Library/Frameworks/R.framework/Versions/4.0/resvources/library")
#remove.packages("Seurat", lib="/Library/Frameworks/R.framework/Versions/4.0/resvources/library")


plotMA( resv, ylim = c(-1, 1) )
plotDispEsts( ddsv )
hist( resv$pvalue, breaks=20, col="grey" )

###  throw out lowly expresvsed genes?? ... I leave this as an exercise
###  add external annotation to "gene ids"
# log transform
rldv = rlog( ddsv )
## this is where you could just extract the log transformed normalized data for each sample ID, and then roll your own analysis
head( assay(rldv) )
mydata = assay(rldv)

sampleDistsv = dist( t( assay(rldv) ) )

#install R 3.5 to avoid maybe the seurat and deseq2 conflict
# heat map
sampleDistMatrixv = as.matrix( sampleDistsv )
rownames(sampleDistMatrixv) = rldv$TissueCode
colnames(sampleDistMatrixv) = NULL
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrixv, trace="none", col=colours)

#turn off developer graphics
dev.off()

# PCs
# wow you can sure tell tissue apart
print( plotPCA( rldv, intgroup = "TissueCode") )
# heat map with gene clustering
library( "genefilter" )
# these are the top genes (that tell tissue apart no doubt)
topVarGenesv <- head( order( rowVars( assay(rldv) ), decreasing=TRUE ), 20 )

png(filename='LamHeatMap_AnnotOnly.png', width=800, height=750)
heatmap.2( assay(rldv)[ topVarGenesv, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
# volcano plot this is an exercise
graphics.off()

AnnotatedGenesInterest<-assay(rldv)[ topVarGenesv, ]
write.csv(AnnotatedGenesInterest,"AnnotatedGenesInterest.csv",row.names = TRUE)

png(filename='LKLDVolcanoAnnotOnly.png', width=800, height=750)
EnhancedVolcano(resv,
                lab = rownames(resv),
                x = 'log2FoldChange',
                y = 'pvalue')
graphics.off()



