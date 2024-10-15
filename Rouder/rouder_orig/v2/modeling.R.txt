library(R2jags)


makeJagsDat=function(dat){
  y=dat$y
  task=dat$task
  sub=dat$sub
  N=length(y)
  I=max(sub)
  J=max(task)
  return(list("y"=y,"sub"=sub,"task"=task,"N"=N,"I"=I,"J"=J))}


estM2=function(dat,M=1000){
  startMu=tapply(dat$y,list(dat$sub,dat$task),mean)
  inits=list(list(
    mu=startMu,
    pSig=rep(1/15^2,dat$J),
    pDel=1/20^2,
    nu=75))
  parameters <- c("mu","pSig","nu", "pDel")
  samples <- jags(dat, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "model2.txt", 
                  n.chains=1,n.iter=M,n.burnin=10,n.thin=1)
  return(samples)
}

estM1=function(dat,M=1000){
  startMu=tapply(dat$y,list(dat$sub,dat$task),mean)
  inits=list(list(
    mu=startMu,
    pSig=rep(1/15^2,dat$J),
    pDel=1/20^2,
    nu=75))
  parameters <- c("mu","pSig","nu", "pDel")
  samples <- jags(dat, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "model1.txt", 
                  n.chains=1,n.iter=M,n.burnin=10,n.thin=1)
  return(samples)
}


estM2=function(dat,M=1000){
  startMu=tapply(dat$y,list(dat$sub,dat$task),mean)
  inits=list(list(
    mu=startMu,
    pSig=rep(1/15^2,dat$J),
    pDel=1/20^2,
    nu=rep(75,dat$J)))
  parameters <- c("mu","pSig","nu", "pDel")
  samples <- jags(dat, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "model2.txt", 
                  n.chains=1,n.iter=M,n.burnin=10,n.thin=1)
  return(samples)
}


estM3=function(dat,M=1000){
  startMu=tapply(dat$y,list(dat$sub,dat$task),mean)
  ppl=tapply(dat$y,dat$sub,mean)
  startPhi=(ppl-mean(ppl))/sd(ppl)
  inits=list(list(
    mu=startMu,
    pSig=rep(1/15^2,dat$J),
    pDel=1/20^2,
    nu=rep(75,dat$J),
    lambda=rep(20,dat$J),
    phi=startPhi))
  parameters <- c("mu","pSig","nu","lambda","phi","pDel")
  samples <- jags(dat, 
                  inits=inits, 
                  parameters=parameters, 
                  model.file = "model3.txt", 
                  n.chains=1,n.iter=M,n.burnin=10,n.thin=1)
  return(samples)
}



