# Downloads expression matrix data from GEO, saves it to .csv files
# Change the filepaths below to change the input & output directories
library(GEOquery)
library(SRAdb)
library(Biobase)
library(readxl)

#Change these directories to customize input/output filepaths
output_directory <- "C:/Users/amusc/Dropbox/Blood Project/Microarray Expression Data/"
metadata_filepath <- "C:/Users/amusc/Dropbox/Blood Project/GEO_Metadata.csv"

#Read accession number list from metadata file
metadata <- read.csv(metadata_filepath)
accession_numbers <- as.vector(metadata[,1])

#Iterate through accession numbers, calling Bioconductor download function to download expression data
for (i in 1:(length(accession_numbers))) {
  print(accession_numbers[i])
  gsm <- getGEO(accession_numbers[i])
  data <- as.data.frame(Table(gsm))
  write.csv(data, paste(output_directory, accession_numbers[i], "_expression_matrix.csv", sep=""))
}
