library(data.table)

# Function to import 10X data that is Gzipped
read10X.gz <- function(folder, id = NULL, transpose = FALSE){
  
  if(is.null(id)) id <- basename(folder)
  
  m <- fread(paste0("zcat < ", folder, "/", "matrix.mtx.gz"), header = TRUE)
  # First row is dimension so append those
  header <- as.numeric(colnames(m))
  sm <- sparseMatrix(i = c(m[[1]], header[1]), j = c(m[[2]], header[2]), x = c(m[[3]], 0))
  rm(m)
  if(transpose) sm <- t(sm)
  colnames(sm) <- paste0(fread(paste0(folder, "/", "barcodes.tsv"), header = FALSE)[[1]], "-", id)
  rownames(sm) <- make.unique(fread(paste0(folder, "/", "genes.tsv"), header = FALSE)[[2]])
  return(sm)
}

LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}
