setwd("/path/to/your/directory")
directory <- getwd()
library("DESeq2")
library("RColorBrewer")
library("gplots")

#for this example I used "under 30" and "over 30" in my file names just as a quick and gross way to sort
files <- grep("30", list.files(directory), value=TRUE)
condition <- (c("below30","over30"))

#make a metadata table, load read count fies and run deseq and make a results matrix
Table <- data.frame(sampleName=files,fileNames=files,condition=condition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable=Table,directory=directory,design=~condition)
dds <- DESeq(dds)
res <- results(dds)

#do a variance stabilizing transformation since they come from seperate patients with different mean count values
vsd <- varianceStabilizingTransformation(dds)

#order by count values and select however many genes you want to look at
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
select <- select[1:100]

#gross heatmap with gross color pallette, but it gives the general gist of what you could do to better visualize it
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(vsd)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale="none",dendrogram="none", trace="none")
