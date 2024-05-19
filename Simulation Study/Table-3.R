n.seq      = c(150,300,600,1000)
beta.seq   = c(1,2,3)
lambda.seq = c(2,3)
alpha.seq  = c(1,2)
resultados = c()


for(b in 1:length(beta.seq))
{
  for(c in 1:length(lambda.seq))
  {
    for(d in 1:length(alpha.seq))
    {
      res.n=c()
      for(a in 1:length(n.seq))
      {
        n         = n.seq[a]
        beta      = beta.seq[b]
        lambda    = lambda.seq[c]
        alpha     = alpha.seq[d]
        nombre    = paste("casosbeta",beta,"lambda",lambda,"alpha",alpha,"n",n,".txt",sep="")
        base      = read.table(nombre, h=F)
        real      = c(beta,lambda,alpha)
        sesgo.aux = base[,1:3]-matrix(real, ncol=3, nrow=nrow(base),byrow=T)
        sesgo     = apply(sesgo.aux, 2, mean)
        se        = apply(base[,4:6], 2, mean)
        RMSE      = sqrt(apply(sesgo.aux^2, 2, mean))
        cp.aux    = matrix(0, ncol=3, nrow=nrow(base))
        for(m in 1:nrow(base))
        {
          for(jj in 1:3)
          {
            lim.inf      = sesgo.aux[m,jj]-1.96*base[m,jj+3]
            lim.sup      = sesgo.aux[m,jj]+1.96*base[m,jj+3]
            cp.aux[m,jj] = ifelse(lim.inf<0 & lim.sup>0, 1, 0)
          }
        }
        cp    = apply(cp.aux, 2, mean)
        res.n = cbind(res.n, cbind(sesgo, se, RMSE, cp))
      }
      resultados=rbind(resultados, res.n)
    }	
  }
}


colnames = rep(c("bias","SE","RMSE","CP"),3)
rownames = rep(c("beta","lambda","alpha"),nrow(resultados)/3)

round(resultados, 4)
table(round(resultados, 4))