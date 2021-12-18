# Розв'язати задачу Коші y' = f(x,y); y(x0) = y0
#методами Ейлера та Рунге Кутта
#варіант 08


A <- 1
y0 <- -2
h <- 0.4
I <- c(0,4)
x0 <- I[1]
xn <- I[2]
n <- length(seq(x0,xn,h))
func <- function(x,y,A){A/(0.4*y^2+2*x)}

# Eiler
x <- seq(x0,xn,h)
y<- y0
for (i in 2:n){
  y <- c(y, y[i-1] + h*func(x[i-1],y[i-1],1))
}
# print(y)
Eiler <- data.frame(x,y)


# Runge-Kutta

yr <- y0
k <- matrix(0,11,4)
for (i in 2:n){
  k[i-1 ,1] <- h* func(x[i-1], yr[i-1],1)
  k[i-1 ,2] <- h * func(x[i-1]+h/2,yr[i-1] + k[i-1 ,1]/2,1)
  k[i-1 ,3] <- h * func(x[i-1]+h/2,yr[i-1] + k[i-1 ,2]/2,1)
  k[i-1 ,4] <- h * func(x[i-1]+h,yr[i-1] + k[i-1 ,3],1)

  yr <- c(yr, y[i-1]+sum(k[i-1 ,1], 2*(k[i-1 ,2:3]), k[i-1 ,4])/6)
}

Runge <- data.frame(x,yr)
cat("Eiler\n")
print(Euklid)
cat("\nRunge\n")
print(Runge)
