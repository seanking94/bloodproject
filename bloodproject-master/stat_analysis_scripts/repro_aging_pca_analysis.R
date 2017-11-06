source("http://www.bioconductor.org/biocLite.R")
library(limma)
library(gplots)
library(ggfortify)

#change this to your working directory
setwd("your/working/directory/")

#fill this in with your potential profile
potential_profile <- c("ENST00000617582.1",
                       "ENST00000485816.5",
                       "ENST00000453218.6",
                       "ENST00000639924.1",
                       "ENST00000532575.1",
                       "ENST00000493679.5",
                       "ENST00000301194.8",
                       "ENST00000442580.5",
                       "ENST00000546844.1",
                       "ENST00000379115.9",
                       "ENST00000496661.1",
                       "ENST00000536616.5",
                       "ENST00000523456.1",
                       "ENST00000501272.6",
                       "ENST00000496835.6",
                       "ENST00000558670.1")

#import expression matrix and metadata files
expr_data <- as.data.frame(read.table("v8_RSEMv1.3.0_transcript_tpm_wholeblood_profile_norm.txt", sep = "\t", header=TRUE, check.names = FALSE))
metadata <- read.table("subject_table.txt", header=TRUE)

#filter metadata to only females between 20 & 50
female_only <- metadata[metadata$SEX == 2, ]
female_only_age_20_to_50 <- female_only[female_only$AGE <= 50, ]
female_only_age_20_to_50 <- female_only_age_20_to_50[!is.na(female_only_age_20_to_50$SUBJID),]
rownames(female_only_age_20_to_50) <- female_only_age_20_to_50$SUBJID

#cross-reference metadata and expression matrices
filtered_expr_data <- expr_data[,colnames(expr_data) %in% female_only_age_20_to_50$SUBJID]
female_only_age_20_to_50 <- female_only_age_20_to_50[female_only_age_20_to_50$SUBJID %in% colnames(expr_data), ]

#ensure rows of expression matrix are ordered the same as columns of metadata
female_only_age_20_to_50 <- female_only_age_20_to_50[order(female_only_age_20_to_50$SUBJID),]
filtered_expr_data <- filtered_expr_data[,order(names(filtered_expr_data))]

#filter data so it contains only genes in your potential profile
filtered_expr_data <- filtered_expr_data[rowSums(filtered_expr_data) != 0,]
profile <- filtered_expr_data[rownames(filtered_expr_data) %in% potential_profile,]

#PCA
pca <- prcomp(t(profile), center = TRUE, scale. = TRUE)
labels <- as.data.frame(female_only_age_20_to_50$AGE)

#summary of explained variances
summary <- summary(pca)
exp_variance <- summary$importance

#plot PCA, change color values in scale_color_gradient if you want to
p <- autoplot(pca, data = labels, colour=colnames(labels)[1], main="PCA of GTEx Samples Using 16 Potential Reproductive Aging Profile Genes", xlab="PC1 (20.07% of variance)", ylab="PC2 (9.491% of variance)")
p + scale_colour_gradient(low = "#143eec", high = "#ec2014", space = "Lab", na.value = "grey50", guide = "colourbar")

#find correlation of first principle component with ages of samples
PCs <- as.data.frame(pca$x)
cor(PCs$PC1, female_only_age_20_to_50$AGE)
