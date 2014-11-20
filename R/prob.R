#' Probability of answer given item and theta.
#' 
#' Returns a vector (for G3PLM) or a matrix (GRM,SM or GPCM) of probabilities for one or more items,
#' given a specified ability (vector).
#' 
#' Detailed description, response models.
#' 
#' @param theta Numeric, vector of length Q ability parameters.
#' @param items Itembank, or a subset of items from an itembank. See also \code{\link{genItembank}}
#' @param model Character, IRT response model. Can be any of "G3PL", "GRM", "SM" or "GPCM".
#' @return List with two entries, a vector or matrix \code{p} of probabilites per item/category, 
#' and a vector \code{p0} of probabilities defined as \code{p0 = q = 1 - sum(P_ij)} 
#' 
prob <- function(theta=0,items=genItembank(model=model),model="GPCM"){
  # TODO: Check input.
  
  # logistic function
  lf <- function(x){ exp(x)/(1+exp(x)) }
  
  # set up output matrix
  P <- matrix(NA,nrow(items$beta),ncol(items$beta)+1)
  
  m <- items$m
  a <- items$alpha
  b <- items$beta
  K <- nrow(b)
  #print(c(a,theta))
  at <- apply(a * drop(theta),1,sum)
  
  # Three Paramater Logistic (MultiDim) (Segall, 1996)
  if(model=="G3PLM"){
    c <- items$guessing
    aux <- numeric(K)
    for (i in 1:K){
      aux[i] <- -a[i,] %*% (theta - b[i,])
    }
    P[,2] <- c + (1-c)/(1+exp(aux))
    P[,1] <- 1 - P[,2]
  }
  
  # graded response model (Samejima, 1969)
  if(model=="GRM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      P[i,1:(m[i]+1)] <- Psi[1:(m[i]+1)] - Psi[2:(m[i]+2)]
    }
  }
  
  # Sequential Model (Tutz, 1990)
  if(model=="SM"){
    for(i in 1:K){
      Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
      for(j in 1:(m[i]+1)){
        P[i,j] <- prod(Psi[1:j]) * (1 - Psi[j+1])
      }
    }
  }
  
  # Generalised Partial Credit Model (Muraki, 1992)
  if(model=="GPCM"){
    for(i in 1:K){
      aux <- exp((1:m[i]) * at[i] - b[i,1:m[i]])
      aux2 <- 1 + sum(aux)
      P[i,2:(m[i]+1)] <- aux / aux2
      P[i,1] <- 1-sum(P[i,],na.rm=TRUE)
    }
  }
  return(P)
}


#### Everything below here is testing code. TODO: remove!
testit <- function(model="GRM",alpha=1,extra=F,s.plot=F){
  x <- seq(-3,3,.1)
  item <- genItembank(Q=1,K=1,M=list(min=3,max=3),model=model)
  
  # set up a fixed item.
  item$alpha <- matrix(alpha)
  item$beta <- matrix(c(-2,-2,0),1)
  item$eta <- matrix(c(-2,0,2),1)
  if(model!="GPCM"){
    item$beta <- item$eta
  }
  if(model=="G3PLM") item$beta <- matrix(0)
  
  y  <- y2 <- y3 <- matrix(0,length(x),4)
  
  for (i in 1:length(x)){
    y[i,] <- prob(x[i],item,model=model)
    if(model =="G3PLM" & !s.plot) y2[i,] <- Pi(x[i],matrix(c(item$alpha,item$beta,item$guessing,1),ncol=4))$Pi
    if(all(model!="SM",model!="G3PLM",!s.plot)) y2[i,] <- Pi(x[i],matrix(c(item$alpha,item$eta),ncol=4),model=model)$Pi
    if(model=='GPCM' & !s.plot) y3[i,] <- unlist(p.(x[i],list(a=item$alpha,b=item$beta,m=3)))
  }
  
  par(mfrow=c(switch(model,"GPCM"=3,"SM"=1,2),1))
  xlab = "\theta"
  if(s.plot) par(mfrow=c(1,1))
  matplot(x,y,type='l',main=paste("MCAT",model,sep=" - "))
  if(model!="SM"& !s.plot) matplot(x,y2,type='l',main=paste("catR",model,sep=" - "))
  if(model=='GPCM'& !s.plot) matplot(x,y3,type='l',main=paste("MCAT(old)",model,sep=" - "))
  
  if(extra) print(ex <- cbind(y[,1],y2[,1],y3[,1]))
  if(extra) all(apply(ex,1,function(x) {
    
  }))
}

testit("G3PLM",1,F,T)
# 
# Q <- 1
# model <- "G3PLM"
# items <- genItembank(model=model,Q=Q,K=5)
# prob(theta=rep(0,Q),model=model,items=items)
# item <- genItembank(model="G3PLM",Q=1,K=1)