#Kholetskii method

set_for_manrix <- c(14.11, 4.45,6.97,3.23,6.66,14.37,13.23,4.18,6.43,12.27,13.67,3.88,3.77,6.22,2.86,11.90)
KHMatrix <- print(t(matrix(set_for_manrix,ncol= 4, nrow=4)))
b<- c(-3.19,-9.30,3.83,19.43,8.16,19.45)
n <- 4
L <- matrix(0,nrow=n,ncol=n)
U <- matrix(0,nrow=n,ncol=n)

# L i1 elements
for ( i in 1:n){
  L[i,1] <- KHMatrix[i,1]
}
# U 1j elements
for (j in 1:n){
  U[j,j]<- 1
  U[1,j] <- KHMatrix[1,j]/L[1,1]
}


for (i in 2:n){
  for( j in 2:n){

    s <- 0
    if (i<j){
      for (k in 1:(i-1)){
        s<- s + L[i,k]*U[k,j]
      }
      U[i,j] <- round((KHMatrix[i,j]-s)/L[i,i],3)
    }
    if ( i>= j){
      for( k in 1:(j-1)){
        s <- s + L[i,k]*U[k,j]
      }
      L[i,j] <- round(KHMatrix[i,j] - s,3)
    }
  }
}


cat(' Matrix A:\n')
print(KHMatrix)
cat(' Matrix U:\n')
print(U)
cat(' Matrix L:]\n')
print(L)
print(det(U)*det(L))
print(det(KHMatrix))


# RUN Method

d <- c(6.8,7.5,5.6,12.2)
a <- c(0,3.5,4.6,6.6)
b <- c(11.5,9.9,9.3,8.7)
c <- c(-2.8,-3.5,-1.7,0)
y <- c(b[1])



alp <- c(-c[1]/y[1])
bet <- c(d[1]/y[1])

for (i in 2:(length(d))){
  y <- c(y, b[i] + a[i]*alp[i-1])
  alp <-c(alp, -c[i]/y[i])
  bet <- c(bet, (d[i]-a[i]*bet[i-1])/y[i])
}
cat('Ps',alp)
cat('Qs', bet)

# y<- c(y, b[length(d)]+a[length(d)]*alp[length(d)-1])
# bet <- c(bet, (d[length(d)]-a[length(d)]*bet[length(d)-1])/y[length(d)])
x<-c(bet[length(d)])


for (i in (length(d)-1):1){
  x<- c(alp[i]*x[1]+bet[i],x)
}
cat( 'Run Method\n')
cat('A,B and Cs:\n', a,'\n',b,'\n',c)
cat('\nD: \n', d )
cat( '\n Alphas: \n', alp)
cat('\nBettas: \n', bet)
cat('\nY\n', y)
cat('\n X: \n', x)
