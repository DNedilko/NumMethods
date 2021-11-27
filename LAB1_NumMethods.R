
library(stringr)
library(Deriv)

func <- function (x){
  return (0.5^(2*x+1))
}

concL <- function (i,j, xs){
  if (i!=j){
    return (paste("(( x - ",xs[j],")/(", xs[i],"-",xs[j],"))"))
  }
  else{
    return (1)
  }
}

funsL <- function (x,i,j,xs){
  if (i!=j){
    return ((( x - xs[j])/(xs[i] - xs[j])))
  }
  else{
    return (1)
  }
}

Lagr <- function(x,xs=seq(1,5,by=1)){
  return (sum(sapply(1:n, function(i) ys[i]*(prod(sapply(1:n, function(j)  funsL(x,i,j,xs)))))))
}

# q <- function (x){
#   return ((x-1)/0.5)
# }

concN <- function(n,y){
  seq <- paste(sapply(2:n, function (i) paste0("(t-",i-2,")",sep="")),collapse = "*")
  polin <- paste(seq, paste( y, factorial(n-1), sep = "/"),sep = "*")
  return (polin)
}

funsN <- function(n,y,x){
  t <- x-1
  vals <- prod(sapply(1:(n-1), function(i) (t-(i-1))))*y/factorial(n-1)
  return (vals)

}

# Newt <- function(ys,)

n<-5
xs <- seq(1,5,by=1)
ys <- sapply(xs,func)
# Значення в точках многочлена Лагранджа
lgj <- sapply(xs, function (x) Lagr(x))

#аналітичний вираз Лагранджа
lagan <- paste(sapply(1:n, function(i) paste(ys[i],paste((sapply(1:n, function(j) concL(i,j,xs))),collapse = " * "),sep=" * ")), collapse="\n + ")
cat("Lagrange polynom \n Ln(x) = ", lagan)
tosimp <- Simplify(lagan, env = parent.frame(), scache = new.env())
cat("Simplified:\n Ln(x) = ",tosimp, "\n Getting rid of brakets: \n Ln(x) = x^4/192 - 7x^3/96 + 71x^2/192 - 77x/96 + 5/8")

#Многочлен Ньютона
delys <- sapply(1:4,function (i) ys[i+1]-ys[i])
del2ys <- sapply(1:3,function (i) delys[i+1]-delys[i])
del3ys <- sapply(1:2,function (i) del2ys[i+1]-del2ys[i])
del4ys <- del3ys[2]-del3ys[1]

fornewys <- c(ys[1],delys[1],del2ys[1],del3ys[1],del4ys)

i <- 10
hn <- 0.5
xn <- seq(1,5.5,hn)

# Аналітичний вираз многочлена Ньютона
newan <- paste(fornewys[1],paste(sapply(2:n, function(x) concN(x,fornewys[x])),collapse = "+"),sep = "+")
cat("\n ---------\nNewton polynom\n N(t) = ", newan, '\n Where t is x-x0, x0=1\n------\nN(x) = ', str_replace_all(newan,"t","x-1" ))

tosimp2 <- Simplify(str_replace_all(newan,"t","x-1" ), env = parent.frame(), scache = new.env())
cat("\nSimplified:\n",tosimp2, "\nGetting rid of brakets: \n 0.38134765625 + 0.00164794921875*x^4 + 0.14556884765625*x^2 - 0.0252685546875*x^3 - 0.3782958984375*x")

# Значення многочлена Ньютона функції у точках
newval <- sapply(xn,function(val) round(sum(fornewys[1],sapply(2:n, function(i) funsN(i,fornewys[i], val))),10))
valsy <- sapply(xn, func)
lagval <- sapply(xn, function (x) Lagr(x))

nwt <- sapply(xs,function(val) round(sum(fornewys[1],sapply(2:n, function(i) funsN(i,fornewys[i], val))),10))

cat('\n ------For 10 values-----\n',xn)
cat('\n Value of function \n',valsy)
cat('\n Newton value \n',newval)
cat('\n Lagrange value \n',lagval)

cat('\n ------For 5 values------\n',xs)
cat('\n Value of function \n',ys)
cat('\n Newton value \n',nwt)
cat('\n Lagrange value \n',lgj)

library(ggplot2)

# Building Graphs
cat("\n Building graph for 5 values")
plot(ys,xs,type="o", col="blue", pch="o", )

points(nwt,xs, col=2, pch="*")
lines(nwt,xs, col=2,lty=2)

points(lgj,xs, col=3,pch="+")
lines(lgj,xs, col=3, lty=3)

legend(8, 0.1, legend=c("Function", "Newton", "Lagrange"), col=c("blue",2,3), lty=1:2, cex=0.8)
text(ys+0.002,xs+0.05,xs)
#
# + plot(valsy,type="b", col=2) + plot(valsy,type="b", col=3)
dots_for_6 <- c(3.7,8)
voffd <- sapply(dots_for_6,func)
newfd <- sapply(dots_for_6,function(val) round(sum(fornewys[1],sapply(2:n, function(i) funsN(i,fornewys[i], val))),10))
lagrfd <- sapply(dots_for_6, function (x) Lagr(x))

cat("\n Value of function for choosen dots 3.7 and 8\n", voffd,
    "\n Newton for 3.7 and 8 \n", newfd, "\n Lagrange for 3.7 and 8\n", lagrfd)

library(caracas)

# n <-



# simpleCondition()
# sq <- sapply(1:n, function(i) paste(ys[i],args[i,],sep=" * "))
# print(sq)
# pol <- sapply(ys,function(x) paste(x,))
#
# Lx <- print()