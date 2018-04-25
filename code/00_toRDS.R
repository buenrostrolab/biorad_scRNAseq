library(data.table)
library(R.utils)

BR_to_10X <- function(file, name){
  x <- data.frame(fread(paste0("zcat < ",file)))
  m <- Matrix::Matrix(data.matrix(x[,-1]))
  
  folder <- paste0("../data/10x_gzip/", name)
  dir.create(file.path(folder), showWarnings = FALSE)
  
  # make data.frames
  geneDF <- data.frame(genes1 = colnames(x)[-1], genes2 =colnames(x)[-1])
  sampleDF <- data.frame(samples = x$CellId)
  
  # write files
  Matrix::writeMM(m, file = paste0(file.path(folder), "/", "matrix.mtx"))
  write.table(geneDF, file = paste0(file.path(folder), "/", "genes.tsv"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  write.table(sampleDF, file = paste0(file.path(folder), "/", "barcodes.tsv"),
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  gzip(paste0(folder,  "/", "matrix.mtx"))
  name
}

df <- data.frame(
  file = list.files("../data/raw/", full.names = TRUE), 
  name = gsub("_S1", "", gsub( ".counts.umiCounts.passingKneeFilter.table.csv.gz", "",list.files("../data/raw/"))), 
  stringsAsFactors = FALSE
)

lapply(1:dim(df)[1], function(i){
  file <- df[i, 1]
  name <- df[i,2]
  BR_to_10X(file, name)
})

