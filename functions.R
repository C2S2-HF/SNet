library(BeSS)
library(plyr)
library(igraph)

BlockGen <- function(V,K,inter.p=0.8,between.p=0.2){
    nodes <- rep(1:V)
    membership <- sample(K,size=V,replace=TRUE)
    node1 <- node2 <- NULL
    for(i in 1:(V-1)){
        for(j in (i+1):V){
            if(membership[i]==membership[j]){
                connect <- rbinom(1,1,inter.p)
                if(connect==1){
                    node1 <- c(node1,i)
                    node2 <- c(node2,j)
                }
            }else{
                connect <- rbinom(1,1,between.p)
                if(connect==1){
                    node1 <- c(node1,i)
                    node2 <- c(node2,j)
                }

            }
        } 
    }
    EdgeList <- data.frame(node1=node1,node2=node2)
    list(Edge=EdgeList,Block=membership)
}

AlphaGen.from.Block <- function(Block,s=4,scale=1){
    K <- length(unique(Block))
    V <- length(Block)
    centroids <- seq(-scale,scale,length.out=K)
    alpha <- rep(0,V)
    for(v in 1:V){
        alpha[v] <- rnorm(1,mean=centroids[Block[v]],sd=s)
    }
    shift <- mean(alpha)
    alpha <- alpha - shift
    centroids <- centroids - shift
    list(alpha=alpha,centroids=centroids)
}

Regression.Gen <- function(alpha,p,sigma=1,beta.given=NULL,correlated=NULL){
    if(is.null(correlated)){
      if(is.null(beta.given)){
      beta <- matrix(rnorm(p,mean=1),ncol=1)
      n <- length(alpha)
      sig<-matrix(0.5,p,p)+diag(p)*0.5
      X<-mvrnorm(n=n,rep(0,p),sig)
      EY <- X%*%beta + alpha
      Y <- EY + rnorm(n)*sigma
      }else{
      beta <- matrix(beta.given,ncol=1)
      n <- length(alpha)
      sig<-matrix(0.5,p,p)+diag(p)*0.5
      X<-mvrnorm(n=n,rep(0,p),sig)
      EY <- X%*%beta + alpha
      Y <- EY + rnorm(n)*sigma
      }
    }
    else{
      if(is.null(beta.given)){
      beta <- matrix(rnorm(p,mean=1),ncol=1)
      n <- length(alpha)
      X <- matrix(rnorm(n*p),ncol=p)
      X <- scale(X,center=TRUE,scale=FALSE)
      EY <- X%*%beta + alpha
      Y <- EY + rnorm(n)*sigma
      }else{
      beta <- matrix(beta.given,ncol=1)
      n <- length(alpha)
      X <- matrix(rnorm(n*p,sd=5),ncol=p)
      X <- scale(X,center=TRUE,scale=FALSE)
      EY <- X%*%beta + alpha
      Y <- EY + rnorm(n)*sigma
      }
    }
    list(X=X,Y=Y,alpha=alpha,beta=beta,EY=EY)
}

lm.net <- function(X,Y,Adj,lambda,theta,cv=NULL,normal=FALSE){
    n <- nrow(X)
    p <- ncol(X)
    Y <- matrix(Y,ncol=1)
    D <- diag(rowSums(Adj))
    L <- D - Adj
    Omega <- lambda*(L + theta*diag(rep(1,n)))
    if(normal){
      d.seq <- rowSums(Adj)
      dd <- diag(1/sqrt(d.seq + theta*rep(1,n)))
      Omega <- lambda * dd%*%(L + theta*diag(rep(1,n)))%*%dd
    }
    X.tilde <- cbind(diag(rep(1,n)),X)
    M <- matrix(0,n+p,n+p)
    M[1:n,1:n] <- Omega
    alpha.beta <- solve(M + t(X.tilde)%*%X.tilde/n,t(X.tilde)%*%Y)/n
    alpha <- alpha.beta[1:n]
    beta <- alpha.beta[-(1:n)]
    cv.MSE <- 0
    if(!is.null(cv)){
        K <- cv
        set.seed(500)
        cv.order <- sample(n,size=n)
        cv.index <- 1:n
        cv.index[cv.order] <- 1:K
        for(k in 1:K){
            current.index <- which(cv.index!=k)
            valid.index <- which(cv.index==k)
            s.A <- Adj[current.index,current.index]
            cv.lm.net <- lm.net(X=matrix(X[current.index,],ncol=ncol(X)),Y=Y[current.index],Adj=s.A,lambda=lambda,theta=theta)
            n.valid <- length(valid.index)
            valid.alpha <- rep(0,n.valid)
            for(v in 1:n.valid){
                valid.id <- valid.index[v]
                valid.alpha[v] <- sum(Adj[valid.id,current.index]*cv.lm.net$alpha)/sum(Adj[valid.id,current.index])
            }
            valid.y <- valid.alpha + X[valid.index,]%*%matrix(cv.lm.net$beta,ncol=1)
            cv.MSE <- cv.MSE + mean((valid.y-Y[valid.index])^2)
        }
        cv.MSE <- cv.MSE/K
    }
    list(alpha=alpha,beta=beta,cv.MSE=cv.MSE)
}

BEss <-function(X,Y,p,kk,imax,cv=NULL){
  n <- nrow(X)
  p <- ncol(X)
  Y <- matrix(Y,ncol=1)
  alpha<-rep(0,n)
  for(j in 1:imax){
    alpha1=alpha
    beta=bess.one(X,Y-alpha1,s=kk,family='gaussian')$beta
    alpha<-Y-X%*%beta
    if(sqrt(sum((alpha1-alpha)^2))<0.001){break}
  }
  cv.MSE <- 0
  if(!is.null(cv)){
    K <- cv
    set.seed(500)
    cv.order <- sample(n,size=n)
    cv.index <- 1:n
    cv.index[cv.order] <- 1:K
    for(k in 1:K){
      current.index <- which(cv.index!=k)
      valid.index <- which(cv.index==k)
      cv.BEss <- BEss(X=matrix(X[current.index,],ncol=ncol(X)),Y=Y[current.index],p=p.x,kk=kk,imax=100,cv=NULL)
      n.valid <- length(valid.index)
      valid.alpha <- rep(0,n.valid)
      for(v in 1:n.valid){
        valid.id <- valid.index[v]
        valid.alpha[v] <- mean(cv.BEss$alpha)
      }
      valid.y <- valid.alpha + X[valid.index,]%*%matrix(cv.BEss$beta,ncol=1)
      cv.MSE <- cv.MSE + mean((valid.y-Y[valid.index])^2)
    }
    cv.MSE <- cv.MSE/K
  }
  list(beta=beta,alpha=alpha,cv.MSE=cv.MSE)
}

SNet <- function(X,Y,Adj,p,kk,imax,lambda,gamma,cv=NULL){
  n <- nrow(X)
  p <- ncol(X)
  Y <- matrix(Y,ncol=1)
  D <- diag(rowSums(Adj))
  L <- D - Adj
  Omega <- 2*n*lambda*(L + gamma*diag(n))
  alpha<-rep(0,n)
  beta<-rep(0,n)
  alpha_mse<-c()
  beta_mse<-c()
  for(j in 1:imax){
    alpha1 = alpha
    beta1 = beta
    beta=bess.one(X,Y-alpha1,s=kk,family='gaussian')$beta
    alpha<-solve(diag(n)+Omega,Y-X%*%beta)
    alpha_mse<-c(alpha_mse,sqrt(sum((alpha1-alpha)^2)))
    beta_mse<-c(beta_mse,sqrt(sum((beta1-beta)^2)))
    if(sqrt(sum((alpha1-alpha)^2))<0.001){break}
  }
  cv.MSE <- 0
  if(!is.null(cv)){
    K <- cv
    set.seed(500)
    cv.order <- sample(n,size=n)
    cv.index <- 1:n
    cv.index[cv.order] <- 1:K
    for(k in 1:K){
      current.index <- which(cv.index!=k)
      valid.index <- which(cv.index==k)
      s.A <- Adj[current.index,current.index]
      cv.SNet <- SNet(X=matrix(X[current.index,],ncol=ncol(X)),Y=Y[current.index],Adj=s.A,p=p.x,kk=kk,imax=100,lambda=lambda,gamma=gamma,cv=NULL)
      n.valid <- length(valid.index)
      valid.alpha <- rep(0,n.valid)
      for(v in 1:n.valid){
        valid.id <- valid.index[v]
        valid.alpha[v] <- sum(Adj[valid.id,current.index]*cv.SNet$alpha)/sum(Adj[valid.id,current.index])
      }
      valid.y <- valid.alpha + X[valid.index,]%*%matrix(cv.SNet$beta,ncol=1)
      cv.MSE <- cv.MSE + mean((valid.y-Y[valid.index])^2)
    }
    cv.MSE <- cv.MSE/K
  }
  list(beta=beta,alpha=alpha,cv.MSE=cv.MSE,amse=alpha_mse,bmse=beta_mse)
}

