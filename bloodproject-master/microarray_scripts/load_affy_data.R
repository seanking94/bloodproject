#set your working directory to where your affy CEL files are
setwd("/path/to/your/directory")
library(affy)

#this automatically loads all the CEL files in the directory
data <- ReadAffy()
affy_frame <- rma(data)

#convert the affy ids to entrez ids
library(hgu133plus2.db)
keys <- hgu133plus2ENTREZID
mapped_probes <- mappedkeys(keys)
mapped_list <- as.list(keys[mapped_probes])

#loop through each row name create corresponding logical vector
logvec <- c()
for (i in 1:length(row.names(affy_frame))){
  if(row.names(affy_frame)[i]%in%labels(mapped_list)){
    logvec[i] <- TRUE
  }
  else{
    logvec[i] <- FALSE
  }
}

#filter the data frame based on the logical vector and keep the ones that have an entrez id
entrez_frame <- affy_frame[logvec,]

#assign entrez ids to each row in the frame, it's affy data so there will be duplicates so append instead of replacing the id
for (i in 1:length(row.names(entrez_frame))){
  row.names(entrez_frame)[i] <- paste0(row.names(entrez_frame)[i], "~" , mapped_list[[row.names(entrez_frame)[i]]])
}
