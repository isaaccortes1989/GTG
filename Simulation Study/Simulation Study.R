rm(list=ls(all=TRUE))
require(actuar)
require(pracma)

rGTG<-function(n,beta,lambda,alpha)
{
  U=runif(n)
  beta*(lambda+qgumbel(U*(1-pgumbel(-lambda,0,1))+pgumbel(-lambda,0,1),0,1))^(1/alpha)
}
#======================================================================#
rGTGn<-function(n,beta,lambda,alpha){
  x <- numeric(n)
  for(i in 1:n){
    x[i]<-rGTG(1,beta,lambda,alpha)
  }
  return(x)
}
#======================================================================#
#n      = 1
#beta   = 1
#lambda = 1
#alpha  = 1
#x      = 2

log.verosimilitud = function(theta, x)
{
  beta=theta[1]
  lambda=theta[2]
  alpha=theta[3]
  logf<--sum(log(alpha)+(alpha-1)*log(x)-alpha*log(beta)-log(1-pgumbel(-lambda,0,1))-(x/beta)^(alpha)+lambda-exp(-((x/beta)^(alpha)-lambda)))
  return(logf)
}


simulacion.GTG<- function(n.iteraciones, tam.muestra,beta,lambda,alpha)
{
  beta.gorro = rep(0, n.iteraciones)
  lambda.gorro = rep(0, n.iteraciones)
  alpha.gorro = rep(0, n.iteraciones)
  media.beta = 0; media.lambda = 0; media.alpha = 0
  sd.beta = 0; sd.lambda = 0; sd.alpha = 0
  iteracion = 0
  while(iteracion < n.iteraciones)
  {
    
    w<- rGTG(tam.muestra,beta,lambda,alpha)
    muestra = w 
    resultado <- optim(par = c(beta,lambda,alpha), fn = log.verosimilitud, 
                       gr = NULL, method = c("L-BFGS-B"), 
                       lower = c(0.0001,-Inf,0.0001), upper = c(Inf, Inf,Inf ),
                       hessian = TRUE, x = muestra)
    if (resultado$convergence == 0)
    {
      iteracion = iteracion + 1
      beta.gorro[iteracion] = resultado$par[1]
      lambda.gorro[iteracion] = resultado$par[2]
      alpha.gorro[iteracion] = resultado$par[3]            
    }
  }
  media.beta = mean(beta.gorro); media.lambda = mean(lambda.gorro) ; media.alpha = mean(alpha.gorro)
  sd.beta = sd(beta.gorro); sd.lambda = sd(lambda.gorro); sd.alpha = sd(alpha.gorro)
  cat("beta.gorro = ", media.beta, "sd.beta = ", sd.beta, fill = TRUE)
  cat("lambda.gorro = ", media.lambda, "sd.lambda = ", sd.lambda, fill = TRUE)
  cat("alpha.gorro = ", media.alpha, "sd.alpha = ", sd.alpha, fill = TRUE)
}

#======================================================================#

n.seq=c(150,300,600,1000)
beta.seq=c(1,2,3)
lambda.seq=c(2,3)
alpha.seq=c(1,2)
replicas=1000
for(a in 1:length(n.seq))
{
  for(b in 1:length(beta.seq))
  {
    for(c in 1:length(lambda.seq))
    {
      for(d in 1:length(alpha.seq))
      {
        n=n.seq[a]
        beta=beta.seq[b]
        lambda=lambda.seq[c]
        alpha=alpha.seq[d]
        resultados=c();m=1
        while(m<=replicas)
        {
          y=rGTGn(n,beta,lambda,alpha)
          #aux=optim(par = c(alpha,delta, beta), fn = log.verosimilitud, 
          #	gr = NULL, method = c("Nelder-Mead"), 
          #	lower = c(0.0001,0.0001,-0.99), upper = c(Inf, Inf,Inf ),
          #	hessian = TRUE, x = y)
          aux=try(optim(par = c(beta,lambda,alpha), fn = log.verosimilitud, 
                        x = y),silent=TRUE)
          
          if(!grepl("Error",aux)[1])
          {
            para=c(aux$par[1],aux$par[2],aux$par[3])
            if(aux$conv==0)
            {
              hes=try(solve(hessian(log.verosimilitud, x0=para, x=y)),silent=TRUE)
              if(!grepl("Error",hes)[1])
              {
                if(all(!is.na(hes)))
                {
                  if(min(diag(hes))>0)
                  {
                    
                    resultados=rbind(resultados, c(para, sqrt(diag(hes))))
                    cc<-paste(c("iteration",m,"beta",beta,"lambda",lambda,"alpha",alpha,"n",n),sep="")
                    print(cc)
                    m=m+1
                  }
                }
              }
            }
          }
        }
        nombre=paste("casosbeta",beta,"lambda",lambda,"alpha",alpha,"n",n,".txt",sep="")
        write(t(resultados),ncolumns=ncol(resultados),file=nombre)
      }
    }	
  }
}