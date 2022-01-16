source(functions.R)
require(doParallel)

registerDoParallel(cores=16)
options(digits=20)
p.x <- 20
set.seed(99)
beta.init <- rep(0,p.x)
beta.trnmb <- sample(1:p.x,5)
beta.init[beta.trnmb]<-rnorm(5)
s.seq <- c(seq(0,3,by=0.2))

data <- foreach(ss=1:16,.packages=c('BeSS','igraph')) %dopar%{
 
  s.value <- s.seq[ss]
  PE.MSE.net <- alpha.MSE.net <- beta.MSE.net <- rep(0,100)
  
  for(m in 1:100){
    set.seed(m+100)
    block.g <- BlockGen(300,K=3,inter.p=0.2,between.p=0.02)
    true.alpha <- AlphaGen.from.Block(block.g$Block,s=s.value,scale=5)
    dd <- Regression.Gen(true.alpha$alpha,p=p.x,sigma=1,beta.given=beta.init)
    g <- graph.data.frame(block.g$Edge,directed=FALSE)
    Adj <- get.adjacency(g,sparse=F)
    new.index <- as.numeric(V(g)$name)
    X <- matrix(dd$X[new.index,],ncol=p.x)
    Y <- dd$Y[new.index]
    EY <- dd$EY[new.index]
    trueAlpha <- dd$alpha[new.index]
    lambda.seq <- seq(0.0001,0,length.out=20)
    mse<-array(0,dim=c(length(lambda.seq),10))
    
    i=1
    for(lambda in lambda.seq){
      for(kk in 1:10){
        PDAS1<-PDAS(X=X,Y=Y,Adj=Adj,p=p.x,kk=kk,imax=500,lambda=lambda,gamma=0,cv=10)
        mse[i,kk]=PDAS1$cv.MSE
      }
      i=i+1
    }
  
    opt.pos <- which(mse == min(mse), arr.ind = TRUE) 
    opt.lambda <- lambda.seq[opt.pos[1,1]]
    opt.k <- opt.pos[1,2]
    opt.model <- PDAS(X=X,Y=Y,Adj=Adj,p=p.x,kk=opt.k,imax=500,lambda=opt.lambda,gamma=0,cv=NULL)    
   
    beta.err.net <- opt.model$beta - dd$beta
    alpha.err.net <- trueAlpha - opt.model$alpha
    PE.err.net <- EY - X%*%matrix(opt.model$beta,ncol=1) - opt.model$alpha
    
    alpha.MSE.net[m] <- sum(alpha.err.net^2)/600
    beta.MSE.net[m] <- sum(beta.err.net^2)/(2*p.x)
    PE.MSE.net[m] <- sum(PE.err.net^2)/600
    
    
  }
  tmp <- list(alpha.MSE.net=alpha.MSE.net,beta.MSE.net=beta.MSE.net,PE.MSE.net=PE.MSE.net)
}

save(data,file="data.Rdata")
