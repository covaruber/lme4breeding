simage <- function(data, Var1=NULL, Var2=NULL, ...){
  M <- table(data[,Var1], data[,Var2])
  M2 <- matrix(M, nrow = nrow(M), ncol = ncol(M))
  name1 <- as.character(substitute(list(Var1)))[-1L]
  name2 <- as.character(substitute(list(Var2)))[-1L]
  Matrix::image(as(t(M2), Class = "dgCMatrix"), 
                xlab=name1, ylab=name2, 
                colorkey=TRUE, ...)
}  