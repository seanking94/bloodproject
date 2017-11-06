library(dglm)

input_file =read.table(input_filename,sep="\t", header=FALSE)
test <-data.frame(genotype = input_file$V1, phenotype = input_file$V2)

out <- dglm(phenotype ~ genotype, ~genotype, data=test, family=gaussian)
s <- summary(out)
pValue <- s$dispersion.summary$coefficients[2, 4]