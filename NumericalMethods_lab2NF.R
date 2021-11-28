library(stringr)
library(numDeriv)


func <- function(x) {3^(x-5)}

runM<- function(xs,ys,h){
  cat("\nMaking coeficients a,b c,d for run method")
  cat("\na\n")
  a <- c(0,rep.int(h,11-2),2*h/3)
  print(round(a,6))
  cat("\nb\n")
  b <- c(2,rep.int(2*h*2,11-2),1+4/3)
  print(round(b,6))
  cat("\nc\n")
  c <- c(1,rep.int(h,11-2),0)
  print(round(c,6))
  cat("\nd\n")
  d <- c(3*(ys[2]-ys[1])/(h*h)-3*grad(func,xs[1])/h,
         sapply(3:10,function(i) 3*((ys[i]-ys[i-1])/h -(ys[i-1]-ys[i-2])/h)), -grad(func,xs[11])-2*(ys[10]-ys[9])/h +3*(ys[11]-ys[10])/h)

  print(round(d,6))
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
  cat("\n as\n")
  as <- ys[1:10]
  print(round(as,6))
  der<- D(D(expression(3^(x-5)),'z'),'z')
  z <- xs[11]
  cat("\n bs\n")
  bst <- (ys[11]-ys[10])/h - (2*cs[9]+cs[10])*h/3
  bs<-c()
  for (i in 1:9){
    btemp<- (ys[i+1]-ys[i])/h - (2*cs[i]+cs[i+1])*h/3
    bs <- c(bs,btemp)
  }
  bs<-c(bs,bst)
  print(round(bs,6))
  z<-xs[1]
  cat("\n ds\n")
  ds <- sapply(1:9,function(i) (cs[i+1]-cs[i])/(3*h))
  dn<- ((ys[11]-ys[10])/h-(ys[10]-ys[9])/h-h*(cs[9]+2*(cs[10]))/3-cs[10]*h)/(h*h)
  ds<- c(ds,dn)
  print(round(ds,6))
  Ss<- sapply(1:10, function(i) cubSpline(as[i],bs[i],cs[i],ds[i],xs[i]))
  # SSs <- c("−103.872553*x^3+25.00722*x^2 − 1.343116*x+0.5","27.596946 * x^3 − 14.43363* x^2+0.704506 * x + 0.481564", "−7.599834*x^3+6.684438*x^2−3.120124*x+0.796375","1.85819*x^3−1.827784*x^2−0.766537*x+0.674168","−0.654901*x^3+1.187926*x^2−2.015756*x+0.945323","0.045843*x^3+0.136809*x^2−1.514861*x+0.96749","−0.151413*x^3+0.49187*x^2−1.942637*x+1.238852","0.017507*x^3+0.137139*x^2−1.133735*x+0.909876","−0.390716*x^3+1.116874*x^2−4.241097*x+3.043063","0.193908*x^3−0.461612*x^2−0.959302*x+1.239506")
  cat("\n PLYNOMS \n")
  print(Ss)
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

a<-1
h<-1
b<-a+10*h
cat("[a,b] -> [",a,",",b,"]")
cat("\n xs")
xs<-print(seq(a,b,h))
cat("\n YS")
ys <- print(sapply(xs,func))
cat("\n Run Method\n")
cs <- round(runM(xs,ys,h),6)
cat("\n value of c\n")
print(cs)
cat("\nIntervals to build spline\n")
intervals<- print(sapply(1:10,function(i) c(xs[i],xs[i+1])))
cat("\n finding a b c d for polynoms")
polynoms<- forSpline(cs,ys,h,xs)
# sapply(polynoms,function(p) cat(p,"\n"))
fgx<- seq(a-1,b+1,0.0001)
fgy<- sapply(fgx,func)
coords_for_splineX <-sapply(1:10,function(i) SplinesDotsX(intervals[1,i],intervals[2,i]))
coords_for_splineY <- sapply(1:10,function(i) SplinesDotsY(coords_for_splineX[,i], polynoms[i]))
# coords_for_splineX <- as.vector(t(splineX))
windows()
# plot(0,0,xlim = c(-10,10),ylim = c(-10,10),type = "n")
plot(fgx,fgy, type='l', col='black')
points(xs,ys, col='red')
k<-701
cols<- rainbow(3)

# setsX<- sapply(1:10,function(i) coords_for_splineX[(1+(i-1)*k):(i*k)])
# setsY<- sapply(1:10,function(i) coords_for_splineY[(1+(i-1)*k):(i*k)])
for (i in 1:10){
  lines(coords_for_splineX[,i],coords_for_splineY[,i],col=cols[i%%2+1], type='l',lty = 3)
}
# lines(as.vector(coords_for_splineX),as.vector(coords_for_splineY), type = "l", lty = 3, col="green")

# library(ggplot2)
# ggplot(data = df, aes(x=x, y=val)) + geom_line(aes(colour=variable))
# line()
