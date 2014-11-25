#' EAP estimate of person parameters.
#' 
#' Estimated-a-posteriori estimate of latent ability. Utilizes Gauss-Hermite quadrature to do the required numerical integration. 
#' Can provide better estimates than ML estimation, particularly early on in a test.
#' 
#' Details
#' 
#' @param theta Numeric (vector), initial theta estimate. Defaults to a null vector of length equal to the number of latent dimensions.
#' @param items Itembank, or a subset of an Itembank. See \code{\link{genItembank}}.
#' @param resp Response pattern, see \code{\link{ans}}.
#' @param model String, valid values are "3PL", "GRM", "SM" or "GPCM". Defaults to "GPCM".
#' @param prior Numeric (matrix), prior covariance.
#' @param debug Logical. Should debug information be shown?
#' @param quad Integer, number of quadrature points \bold{per dimension}. Note that the total number of quadrature points is \code{quad^Q}, where Q is the number of latent dimensions.
#' @return estimate A vector of latent ability estimates.
#' @seealso eapSE.
eapEst <- function(items, resp, model, prior, ip = 4, debug = FALSE,...){
  # initialize Gauss quadrature points
  # TODO: This only needs doing once for the whole test, should happen on test init when any estimator is set to EAP.
  GH <- init.gauss(ip,prior)
  
  # Evaluate that integral!
  # TODO: Somehow this ALWAYS returns around .3, regardless of true theta.
  return(eval.gauss(LL,GH,items=items,resp=resp,model=model,debug=debug))
}


# Likelihood
LL <- function(theta,items,resp,model){
  P <- prob(theta,items,model)
  if(model == "3PL"){
    out <- sum(log(P[,1]^(1-resp)) + log(P[,2]^(resp)))
  } else {
    out <- 0
    j <- resp+1
    for(i in 1:nrow(items$beta)){
      out <- out + log(P[i,j[i]])
    }
  }
  return(out)
}


#' Initialize Gauss-Hermite quadrature points and weights.
#' 
#' Uses fastGHQuad to create unidemsional quadrature points. 
#' These are subsequently expanded across the Q-dimensional grid, and shifted to account for the prior covariance matrix.
#' Weights are computed as the product of weights for each dimension at each quadrature point.
#' The mean vector is assumed to be zero, alternative prior means are not supported at this point.
#' 
#' Details
#' 
#' @param ip Integer, number of quadrature points, \bold{per dimension}. Note that the total number of quadrature points is \code{ip^Q}.
#' @param prior Covariance matrix of the prior distribution.
#' @return X Matrix of dimensions \code{ip^Q} by \code{Q} of quadrature points.
#' @return W Vector of length \code{ip^Q} of weights.
#' @return ip Integer number of quadrature points per dimension used after sanity check.
#' @export
#' @seealso Package fastGHQuad, eval.gauss().
## initialize expanded grid of gauss hermite quadrature points (multidimensional).
init.gauss <- function(ip=10,prior){
  
  Q <- ncol(prior)
    
  # require/install fastGHquad (will be package dependancy).
  # TODO: Implement dependency properly.
  if(!require(fastGHQuad)){ install.packages('fastGHQuad'); require(fastGHQuad)}
  
  # get quadrature points, create grid
  x <- fastGHQuad::gaussHermiteData(ip)
  x$x <- x$x * 2^.5
  # (if anyone knows an easy way to assign a single vector x times to x list elements, please tell me. )
  X <- as.matrix(expand.grid(lapply(apply(replicate(Q,x$x),2,list),unlist)))

  
  # calculate weights
  x$w <- x$w / pi^.5
  # same as above, roundabout way to get the combn weights for each combination of quad points
  g <- as.matrix(expand.grid(lapply(apply(replicate(Q,x$w),2,list),unlist)))
  # combined weight is the product of the individual weights
  W <- apply(g,1,prod)
  
  # compute lambda (eigen decomposition covar matrix)
  lambda <- with(eigen(prior), vectors %*% diag(sqrt(values),Q,Q))
  
  # account for correlation
  X. <- t(lambda %*% t(X))
  
  # TODO: account for mu != 0.
  return(list(X=X.,W=W,ip=ip))
}

#' Evaluation of multivariate normal distributed integral
#' 
#' Evaluates a given function with a (built-in) multivariate normal prior distribution by Gauss-Hermite quadrature.
#'
#' The evaluated function is assumed to have a multivariate normal distribution, with a given mean vector and covariance matrix. 
#' The default identity function \code{function(x) 1} reduces to an integral over a multivariate normal distribution with mean vector \code{mu} and covariance matrix \code{Sigma}.
#' 
#' @param FUN (Likelihood) function of the parameters to be estimated. 
#'     Defaults to \code{funtion(x) 1}, in which case only the built-in multivariate normal pdf is evaluated.
#' @param prior Covariance matrix.
#' @param X Matrix of quadrature points, see \code{\link{init.gauss}}. Alternatively, the list of quadrature points and weights produced by \code{\link{init.gauss}}.
#' @param W Vector of weights, or \code{NULL} if provided in \code{X}.
#' @param ... Additional arguments passed on to FUN.
#' @return A vector with the value evaluated integrals.
#' @seealso \code{\link{init.gauss}} for creating quadrature points.
#' @export
#' @examples
#' quadPoints <- init.gauss(Q=3)
#' # expected value of 3-dimensional multivariate normal distribution: N(0,1). 
#' # (Since mean is currently fixed at zero, this is always zero.)
#' (integral <- eval.gauss(Q=3,X=quadPoints))
#' round(integral)
eval.gauss <- function(FUN = function(x) 0,X=NULL,W=NULL,debug=FALSE,...){
  if (is.list(X)){
    W <- X$W
    X <- X$X
  }
  FUN <- match.fun(FUN)
  if (is.null(X) | is.null(W)) stop("Quadrature points and weights are required. See init.gauss.", call.=F)
  
  ipq <- length(W)
  f <- numeric(ipq)
  
  # main loop
  for (i in 1:ipq){
    f[i] <- FUN(X[i,],...) + log(W[i])
  }
  
  # some numerical safeguards
  m <- 700 - max(f)
  f <- f + m
  
  # back to normal scale
  f <- exp(f)
  
  # normalizing constant
  p1 <- sum(f)
  
  # multiply integrals with x values, sum over columns, divide by normalizing constant.
  out <- colSums(f * X) / p1
  
  # debug plot
  if(debug & ncol(X) == 2){
    ip <- ipq^(1/ncol(X))
    nL <- matrix(f,ip)
    nL <- t(nL[,ncol(nL):1])
    x <- sort(unique(X[,1]))
    y <- sort(unique(X[,2]))
    r <- 2
    contour(x=x,y=y,z=nL,xlim=out[1]+c(-r,r),ylim=out[2]+c(-r,r),axes=FALSE)
    axis(side = 1, at = round(x,2))
    axis(side = 2, at = round(y,2))
    points(out[1],out[2],lwd=2,col='red')
    box()
  }
  if(debug & ncol(X) == 1){
    plot(X,f,type='l')
    print(f)
  }
  
  # return computed integral.
  return(out)
}
