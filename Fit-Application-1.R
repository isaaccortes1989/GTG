rm(list=ls(all=TRUE))

library(tpn)
library(agricolae)
library(fitdistrplus)
library(moments)


z <- read.table("C:/Users/Isaac/Documents/GTG/data.txt",header=F)
x <- z$V11

#============================#
#====Descriptive Analysis====#
#============================#

hist(x)
var(x)
skewness(x)
kurtosis(x)
summary(x)


log.verosimilitud = function(theta, x)
{ 
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  
  beta   = theta[1]
  lambda = theta[2]
  alpha  = theta[3]
  logf   <- -sum(log(alpha)+(alpha-1)*log(x)-alpha*log(beta)-log(1-pgumbel(-lambda,0,1))-(x/beta)^(alpha)+lambda-exp(-((x/beta)^(alpha)-lambda)))
  return(logf)
}


w <- optim(c(1,1,1),log.verosimilitud, x=x, hessian=TRUE, method="L-BFGS-B",lower=c(-Inf,0,-Inf),upper=c(Inf,Inf,Inf))
w$par
sqrt(diag(solve(w$hessian)))
w$value

#======================#
#======Criteria========#
#======================#

AIC = 2*3+2*w$value
BIC = log(length(x))*2+2*w$value
AIC
BIC

EVWeibull<- fitdistr(x=x, densfun='weibull')
EVWeibull
distribucionW <- fitdist(x, distr = 'weibull')
summary(distribucionW)
est.stpn(x)

#======================================#
#ML ESTIMATES==========================#
#======================================#

log.verosimilitud = function(theta, x)
{ 
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  
  beta=theta[1]
  lambda=theta[2]
  alpha=1
  logf<- -sum(log(alpha)+(alpha-1)*log(x)-alpha*log(beta)-log(1-pgumbel(-lambda,0,1))-(x/beta)^(alpha)+lambda-exp(-((x/beta)^(alpha)-lambda)))
  return(logf)
}
w2<-optim(c(1,1),log.verosimilitud,x=x, hessian=TRUE, method="L-BFGS-B",lower=c(0,-Inf),upper=c(Inf,Inf))
w2$par
sqrt(diag(solve(w2$hessian)))

w2$value

AIC2 = 2*2+2*w2$value
BIC2 = log(length(x))*2+2*w2$value
AIC2
BIC2


#=====================#
#=========FIT GTG=====#
#=====================#

dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
qgumbel <- function(p, a, b) a-b*log(-log(p))

Variable <- x
hist(Variable, probability = T, nclass = 10, main = "", xlab = "", ylim=c(0,0.0004))

dggumbel = function(x,sigma,lambda,alpha) ((alpha*x^(alpha-1))/((sigma^alpha)*(1-pgumbel(-lambda,0,1))))*dgumbel(((x/sigma)^alpha-lambda),0,1)

GG=function(sigma,lambda,alpha,a,b,c,d,e,f,g)
{
  x=seq(0.001,10000,0.1)
  density=((alpha*x^(alpha-1))/((sigma^alpha)*(1-pgumbel(-lambda,0,1))))*dgumbel(((x/sigma)^alpha-lambda),0,1)
  density2=((1*x^(1-1))/((a^1)*(1-pgumbel(-b,0,1))))*dgumbel(((x/a)^1-b),0,1)
  density3=dweibull(x,c,d)
  density4=dstpn(x,e,f,g)
  aa=c(expression(GTG),expression(TG),expression(WEI),expression(STPN))
  
  lines(x,density,type="l",col="1")
  lines(x,density2,type="l",col="2")
  lines(x,density3,type="l",col="3")
  lines(x,density4,type="l",col="4")
  legend("topright",aa,lty=c(1,1,1,1),lwd=2,bty="n",cex=2,col=c(1,2,3,4))
}

GG(w$par[1],w$par[2],w$par[3],w2$par[1],w2$par[2],1.607485 ,2454.960160,660.566543,1.973219,2.460686)

