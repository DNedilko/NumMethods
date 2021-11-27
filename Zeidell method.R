seidel <- function(a,b,init,eps=0.00003) {
  n <- ncol(a)
# define initial value for algorithm if it is not supplied;
# in this case, we will define init[i] <- b[i]/a[i,i]
  if (missing(init)) init <- b/diag(a)
  x <- init
  error <- abs(b - a%*%x) # %*% is matrix multiply
  maxerror <- max(error)# iterate until maximum error is smaller than specified threshold eps
  while (max(error)>eps) {
    for (j in 1:n) {
      aj <- as.vector(a[j,]) # j-th row of the matrix
      ax[j] <- (b[j]-sum(aj[-j]*x[-j]))/aj[j]}
    error <- abs(b - a%*%x)
    maxerror <- c(maxerror,max(error))}
  r <- list(x=x,maxerror=maxerror)
  r}

a<-c(20.9,1.2,2.1,0.9,1.2,21.2,1.5,2.5,2.1,1.5,19.8,1.3,0.9,2.5,1.3,32.1)
b<-c(21.7,27.46,28.76,49.72)
seidel(a,b,0)