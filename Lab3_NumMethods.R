# windows()
library(tidyverse)
library(Deriv)
interval <- c(-2.1, 1.1)
a <- interval[1]
b <- interval[2]
func <- function(x) { exp(x)/(2*x+7) }
# plot.function(func)
# Обчислити заданий інтеграл Рімана із точністю e = 10^-5
# ---------------------------------------------------
# 1 - методом трапецій із визначенням
# кількості інтервалів розбиття
# через оцінку похибки

# Nta похідна
Der <- function(fun,x,n){
  fun <- parse(text=fun)
  Def <- D(fun, 'x')
  for (i in 2:n){
    Def <-D(Def,'x')
  }
  ret <- eval(Def, envir=list(x=x))
  return (ret)

}

DerT <- function(fun,t,n){
  fun <- parse(text=fun)
  Def <- D(fun, 't')
  for (i in 2:n){
    Def <-D(Def,'t')
  }
  ret <- eval(Def, envir=list(x=t))
  return (ret)

}

# Value of function in xes
valT<- sapply(interval, function(x) Der('exp(x)/(2*x+7)', x, 2))
cat("\nMaxes = " , valT)
M2 <- max(abs(valT))
cat("\nM2 = ", M2)
eps <- 10^(-5)
h <- sqrt(12*eps/(M2*abs(b-a)))
cat("\nh = ", h)
n <- round(print((b-a)/h) + 0.5)
cat("\nn = ", n)
h <- (b-a)/n
cat( "\nh = ",h)
he <- seq(a, b, length.out = n)
xsT <- func(seq(a, b, length.out = n))
Integral1 <- h*(sum(xsT[1]/2,xsT[n]/2,xsT[2:(n-1)]))
cat( "\ni  |   x   |  f(x)")
for (i in 1:7){ cat(sprintf("\n %d | %.8f | %.8f", i,he[i], xsT[i] ))}
for (i in (length(xsT)-7):length(xsT)) {cat(sprintf("\n %d | %.8f | %.8f", i,he[i], xsT[i] ))}

cat("\n I = ", Integral1)
# ---------------------------------------------------
# 2 - методом Сімпсона із використанням принципу Рунге
q<- 4
k<- 0
check <- TRUE
cat ("\n k | h        | I h        | h/2        | I h/2")
while(check){
  # print(k)
  dil <- 2^k
  ht1 <- h/(dil)
  ht2 <- h/(dil*2)

  nt1 <- round((b-a)/ht1 + 0.5)
  nt2 <- round((b-a)/ht2 + 0.5)

  mt1 <- nt1/2
  mt2 <- nt2/2

  xsSt1 <- seq(a, b, length.out = 2*mt1)
  xsSt2 <- seq(a, b, length.out = 2*mt2)

  S1t1 <- sapply(2:(mt1), function(i) func(xsSt1[2*i-1]))
  S2t1 <- sapply(2:(mt1), function(i) func(xsSt1[2*i-2]))

  S1t2 <- sapply(2:(mt2), function(i) func(xsSt2[2*i-1]))
  S2t2<- sapply(2:(mt2), function(i) func(xsSt2[2*i-2]))

  Integral2t1 <- ((b-a)/(6*mt1))*sum(func(c(a,b)), 4*S1t1,2*S2t1)
  Integral2t2 <- ((b-a)/(6*mt2))*sum(func(c(a,b)), 4*S1t2,2*S2t2)
  epsF <- abs(Integral2t1-Integral2t2)/(2^q - 1)
  cat(sprintf("\n %d | %.8f | %.8f | %.8f | %.8f", k, ht1, Integral2t1,ht2, Integral2t2))
  k <- k+1

  # print(epsF<eps)
  check <- !(epsF<eps)


}

cat("\n I = ", Integral2t2)
# print()
# ---------------------------------------------------
# 3 - за допомогою формул Гаусса
# Поліном Легранджа для моєї функції
n <- 5
v <-c((a+b)/2, (b-a)/2)
xp <- paste(v[1], " + ", v[2],"* t")
funwitht <- str_replace_all( 'exp( x )/(2*( x )+7)',' x ', xp)
dx <- paste(D(parse(text=xp),'t'))
newab <- c(-1,1)
funFI <- paste(dx,paste("(",funwitht,")"),sep=" * ")
cat("\n", xp)
cat("\n", funFI)

ts <- c(-0.90617985, -0.53846931, 0, 0.53846931, 0.90617985)
c <- c(0.23692688, 0.47862868, 0.56888889, 0.47862868, 0.23692688)
ft <-sapply(ts, function (t) eval(parse(text=funFI), envir=list(x=t)))

cat ("\n| t        | f(t)    | c")
for (i in 1:5){
  cat(sprintf("\n| %.5f | %.5f | %.5f",ts[i], ft[i], c[i]))
}

L <- sum(c*ft)
cat("\nI = ", L)


# M10 <- print(sapply(newab, function(t) DerT(funFI,t,10)))
# M10max <- max(M10)
# Rf <- print((M10max*(2^11)*(factorial(5)^4))/((11)*(factorial(10)^3)))
# valst <- sapply(newab, function (x) M10)
# M10mx <- max(M10)

# cat("\n",funFI)


