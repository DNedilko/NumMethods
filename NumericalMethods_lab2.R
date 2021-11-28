library(stringr)
library(numDeriv)


func <- function(x) {0.5^(2*x+1)}

runM<- function(xs,ys,h){
  a <- c(0,rep.int(h,11-2),2*h/3)
  # cat("\nb")
  b <- c(2,rep.int(2*h*2,11-2),1+4/3)
  # cat("\nc")
  c <- c(1,rep.int(h,11-2),0)
  # cat("\nd")
  d <- c(3*(ys[2]-ys[1])/(h*h)-3*grad(func,xs[1])/h,
         sapply(3:10,function(i) 3*((ys[i]-ys[i-1])/h -(ys[i-1]-ys[i-2])/h)), -grad(func,xs[11])-2*(ys[10]-ys[9])/h +3*(ys[11]-ys[10])/h)
  # cat("\n", length(d))
  y <- c(b[1])

  alp <- c(-c[1]/y[1])
  bet <- c(d[1]/y[1])

  for (i in 2:(length(d))){
    y <- c(y, b[i] + a[i]*alp[i-1])
    alp <-c(alp, -c[i]/y[i])
    bet <- c(bet, (d[i]-a[i]*bet[i-1])/y[i])
  }
  x<-c(bet[length(d)])

  for (i in (length(d)-1):1){
    x<- c(alp[i]*x[1]+bet[i],x)
  }

  return (x)
}

cubSpline <- function(a,b,c,d,xi){
  S<- paste(a,"+",b,"*(x - ",xi,")+",c,"*(x - ", xi, ")^2 + ",d,"*(x - ", xi, ")^3")
  return (S)
}



forSpline <- function (cs,ys,h,xs){
  # cat("\n as")
  as <- ys[1:10]
  # print(round(as,3))
  der<- D(D(expression(0.5^(2*z+1)),'z'),'z')
  z <- xs[11]
  # cat("\n bs")
  bst <- (ys[11]-ys[10])/h - (2*cs[9]+cs[10])*h/3
  bs<-c()
  for (i in 1:9){
    btemp<- (ys[i+1]-ys[i])/h - (2*cs[i]+cs[i+1])*h/3
    bs <- c(bs,btemp)
  }
  bs<-c(bs,bst)
  # print(round(bs,4))
  z<-xs[1]
  # cat("\n ds")
  ds <- sapply(1:9,function(i) (cs[i+1]-cs[i])/(3*h))
  dn<- ((ys[11]-ys[10])/h-(ys[10]-ys[9])/h-h*(cs[9]+2*(cs[10]))/3-cs[10]*h)/(h*h)
  ds<- c(ds,dn)
  # print(round(ds,4))
  Ss<- sapply(1:10, function(i) cubSpline(as[i],bs[i],cs[i],ds[i],xs[i]))
  #   # print(Ss)
  return (Ss)
}

SplinesDotsX <- function(interval1,interval2){

  dotsX <- seq(interval1,interval2,0.001)
  return (dotsX)
}
SplinesDotsY<- function(dots, polynom){
  dotsY <- sapply(dots,function(x) eval(parse(text=polynom)))
  return (dotsY)
}

a<-0
h<-0.7
b<-a+10*h
cat("\n xs")
xs<-print(seq(a,b,h))
cat("\n YS")
ys <- print(sapply(xs,func))
cat("\n ew\n")
cs <- print(round(runM(xs,ys,h),4))

cat("\n")
intervals<- sapply(1:10,function(i) c(xs[i],xs[i+1]))

polynoms<- forSpline(cs,ys,h,xs)
# sapply(polynoms,function(p) cat(p,"\n"))
fgx<- seq(a-0.1,b+0.1,0.001)
fgy<- sapply(fgx,func)
coords_for_splineX <-sapply(1:10,function(i) SplinesDotsX(intervals[1,i],intervals[2,i]))
coords_for_splineY <- sapply(1:10,function(i) SplinesDotsY(coords_for_splineX[,i], polynoms[i]))
# coords_for_splineX <- as.vector(t(splineX))
windows()
# plot(0,0,xlim = c(-10,10),ylim = c(-10,10),type = "n")
plot(fgx,fgy, type='l', col='black')
points(xs,ys, col='red')
k<-701
cols<- rainbow(10)

# setsX<- sapply(1:10,function(i) coords_for_splineX[(1+(i-1)*k):(i*k)])
# setsY<- sapply(1:10,function(i) coords_for_splineY[(1+(i-1)*k):(i*k)])
for (i in 1:10){
  lines(coords_for_splineX[,i],coords_for_splineY[,i],col=cols[i], type='l',lty = 3)
}
# lines(as.vector(coords_for_splineX)[1:k],as.vector(coords_for_splineY)[1:k], type = "l", lty = 3, col=c(1,2,3))