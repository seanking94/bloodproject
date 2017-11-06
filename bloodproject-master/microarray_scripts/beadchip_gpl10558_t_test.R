source("http://www.bioconductor.org/biocLite.R")
biocLite(c("beadarray", "limma", "GEOquery", "illuminaHumanv1.db", "illuminaHumanv2.db", 
           "illuminaHumanv3.db", "illuminaHumanv4.db", "BeadArrayUseCases", "GOstats", "GenomicRanges", "Biostrings"))
library(limma)
library(illuminaHumanv4.db)
library(gplots)

#Replace these file paths with the correct ones on the computer you are using
metadata_filepath <- "C:/Users/amusc/Dropbox/Blood Project/GEO_Metadata.csv"
beadchip_expr_filepath <- "C:/Users/amusc/Dropbox/Blood Project/Microarray Expression Data/beadchip_consolidated_expression_data.csv"

#Read expression file and metadata file
metadata <- read.csv(metadata_filepath)
beadchip_expr_data <- read.csv(beadchip_expr_filepath)
rownames(beadchip_expr_data) <- beadchip_expr_data[,1]
beadchip_expr_data <- beadchip_expr_data[,-c(1)]

#filter out non-beadchip expression data
beadchip_metadata <- metadata[metadata$Platform.Assay %in% c("GPL10558", "GPL6947"),]

#had to exclude from analysis due to batch effects
beadchip_metadata <- beadchip_metadata[beadchip_metadata$GSE.. == "GSE58137",]

beadchip_acc_no <- beadchip_metadata$GSM..
rownames(beadchip_metadata) <- beadchip_metadata$GSM..
colnames(beadchip_expr_data) <- beadchip_metadata$GSM..
beadchip_expr_data <- beadchip_expr_data[,colnames(beadchip_expr_data) %in% beadchip_metadata$GSM..]

#separate data into age groups
age_range <- beadchip_metadata[beadchip_metadata$Age <= 25,]
age_range <- rbind(age_range, beadchip_metadata[beadchip_metadata$Age >= 45,])

beadchip_expr_data <- beadchip_expr_data[,colnames(beadchip_expr_data) %in% age_range$GSM..] 

ages <- floor(beadchip_metadata$Age)

group1 <- beadchip_metadata[beadchip_metadata$Age <= 25,]
group1_gsms <- group1$GSM..
group2 <- beadchip_metadata[beadchip_metadata$Age >= 45,]
group2_gsms <- group2$GSM..

ages <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(beadchip_expr_data)))
colnames(ages) <- colnames(beadchip_expr_data)

#label two groups as binary data 
# 0 - age 20-25
# 1 - age 45-50
for (i in 1:ncol(beadchip_expr_data)) {
  id = colnames(beadchip_expr_data)[i]
  if (id %in% group1_gsms) {
    ages[1,i] = 0
  } else {
    ages[1,i] = 1
  }
}


#perform two-sample t-test on expression data, using age as groups
t_test <- as.data.frame(matrix(NA, nrow=nrow(beadchip_expr_data), ncol=1))
rownames(t_test) <- rownames(beadchip_expr_data)

for (i in 1:nrow(beadchip_expr_data)) {
  curr_gene <- beadchip_expr_data[i,]
  curr_gene <- t(curr_gene)
  gr1 <- curr_gene[,colnames(curr_gene) %in% group1_gsms]
  gr2 <- curr_gene[,colnames(curr_gene) %in% group2_gsms]
  t <- t.test(gr1, gr2, paired=FALSE, var.equal = TRUE)
  t_test[i,1] = t
}

#set significance thresholds based on t distribution with n - 2 = (11 + 13) - 2 = 22 degrees of freedom
t_test$'95% Significance' = t_test$V1 >= 2.074
t_test$'99% Significance' = t_test$V1 >= 2.819
t_test$'99.9% Significance' = t_test$V1 >= 3.792

#Map illumina probe IDs to entrez gene IDs and Gene Symbols using limma
symbols <- mapIds(illuminaHumanv4.db, rownames(t_test), "SYMBOL","PROBEID")
entrez <- mapIds(illuminaHumanv4.db, rownames(t_test), "ENTREZID","PROBEID")
t_test[,"Entrez ID"] <- entrez

#filter out unmapped entrez ID's and those genes below t-test significance threshold of 99%
not_na <- t_test[!is.na(t_test$`Entrez ID`),]
not_na <- t_test[!(t_test$`Entrez ID` == "NA"),]
not_na <- t_test[!(t_test$`Entrez ID` == ""),]
not_na <- t_test[t_test$`99% Significance` == TRUE,]
matr = as.matrix(beadchip_expr_data[rownames(beadchip_expr_data) %in% rownames(not_na),])
rownames(matr) <- not_na$`Entrez ID`

#plot significant genes on heatmap, forcing order of samples to be by age group
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
heatmap.2(matr, Colv=FALSE, col=my_palette, dendrogram="row",trace="none", margin=c(8,9), hclust=hclustfunc,distfun=distfunc,symm=F,symkey=F,symbreaks=T, scale="none")
