#' Derivatives
#' 
#' First and second order derivatives of the specified model, for the specified parameters.
#' 
#' Details
#' 
#' @param theta Person parameter.
#' @param items Itembank, see \code{\link{genItembank}}.
#' @param resp Response pattern for the Itembank. Note that all items in the itembank need a corresponding response.
#' @param model IRT model, see \code{\link{MCAT}}. \strong{Defaults to GPCM}.
#' @param prior Covariance matrix if a prior should be included, \code{NULL} otherwise. \strong{Defaults to \code{NULL}}.
#' # Do we have to do something here for information? I.e. don't sum it?
#' @return If \code{sum == TRUE} a list containing a matrix of first order derivatives, 
#' and a list of matrices holding the second order derivatives per item.
#' otherwise a vector with the sum of first order derivatives and a matrix with the sum of second order derivatives. TODO: <- THIS
deriv <- function(theta,items,resp,model="GPCM",prior=NULL){
  a <- items$alpha
  b <- items$beta
  u <- resp
  m <- items$m
  K <- nrow(a)
  Q <- ncol(a)
  d <- D <- numeric(K)
  
  # Set up the various derivatives.
  if (model == "G3PLM"){
    # Segall (MCAT book, 1996, p.71)
    P <- prob(theta,items,model)
    c <- items$guessing
    q <- P[,1]
    p <- P[,2]
    
    d <- ((p-c) * (u-p)) / ((1-c) * p)
    D <- (q *(p-c)*(c*u-p^2)) / (p^2*(1-c)^2)
  }
  
  if (model == "GRM" | model == "SM"){
    # TODO: See if we should put the whole d1, d2 thing in prob. We're recalculating the building blocks here...
    lf <- function(x){ exp(x)/(1+exp(x)) } # from prob.
    at <- apply(a * drop(theta),1,sum)
    
    if (model == "GRM"){
      # Graded Response Model (Glas & Dagohoy, 2007)
      for(i in 1:K){
        j <- u[i]+1 # no more messing about with 0 based indices.
        Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
        d[i] <- 1 - Psi[j] - Psi[j+1]
        D[i] <- -(Psi[j] * (1-Psi[j]) + Psi[j+1] * (1-Psi[j+1]))
      }
    } else {
      # Sequential Model (Tutz, xxxx)
      # TODO: Wording in Glas & Dagohoy for dij (SM) is odd. So Psi_0 = 1, right? Why 1-Psi_(m+1)=1 and not Psi_(m+1)=0?
      # TODO: also, why is there an extra set of parenthesis in the formula for dij?
      # TODO: Finally, it doesn't work. See comments in estiamte.R.
      for(i in 1:K){
        j <- u[i]+1 # no more messing about with 0 based indices.
        Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
        d[i] <- sum((1-Psi[2:j]) - Psi[(2:j)+1])
        D[i] <- -sum(Psi[2:(j+1)] * (1 - Psi[2:(j+1)]))
      }
    }
  }
  
  if (model == "GPCM"){
    # Generalized Partial Credit Model (Muraki, 1992)
    # TODO: again, why the extra parentheses? 
    P <- prob(theta,items,model)
    for(i in 1:K){
      mi <- 1:m[i]
      pi <- P[i,mi+1] # remove j = 0, index is now also correct.
      mp <- sum(mi*pi)
      
      d[i] <- u[i] - mp
      D[i] <- -sum((mi * pi) * (mi - mp))
    }
  }
  
  # Build the derivatives, sum over items.
  # TODO: figure out wether individual derivatives are required anywhere, possibly for information on next item?
  # d1
  res1 <- apply(a * d,2,sum)
  
  # d2
  res2 <- matrix(0,Q,Q)
  for (i in 1:K){
    res2 <- res2 + a[i,] %*% t(a[i,]) * D[i]
  }
  
  # prior
  if (!is.null(prior)) res1 <- res1 - solve(prior)%*%theta
  if (!is.null(prior)) res2 <- res2 - solve(prior)
  
  
  return(list(d1=res1,d2=res2))
}


######## test>
# Q <- 2
# theta <- rep(0,Q)
# model <- "GPCM"
# items <- genItembank(model = model,Q=Q,K=5)
# resp <- ans(theta,items,model)
# solve(deriv(theta,items,resp,model)$d2)
