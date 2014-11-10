#' Generate itembank.
#' 
#' Generates a Q-dimensional itembank under the Graded Response Model (GRM), Sequential Model (SM), Generalized Partial Credit Model (GPCM) or Generalized Three Parameter Logistic Model (G3PLM).
#' 
#' Alpha parameters may be correlated, and are \code{> 0}. 
#' This is accomplished by an initial draw from a multivariate normal distribution with the given covariance matrix. 
#' We then take the exponent of these values, and for each item set the sum of exponents to 1. 
#' \emph{This does mean that the observed covariance is not the same as the generating covariance matrix.}
#' 
#' @param Q integer; number of dimensions
#' @param K integer; number of items
#' @param model character; IRT model used to generate item parameters (can be either 'GRM', 'SM', or 'GPCM").
#' @param covar matrix; Covariance matrix for alpha (discrimination parameter), defaults to identity matrix of rank Q. See \bold{details}.
#' @param a list of options for sampling discrimination parameters. See \bold{details}.
#' @param b list of options for sampling difficulty parameters. See \bold{details}.
#' @param c list(min, max); named list of distribution parameters for the guessing parameter, c.
#'     Guessing parameters are drawn from a continuous uniform distribution with minimum \code{a} and maximum \code{b}. 
#'     Default is list(min=0,max=0) for all \code{c} = 0.
#'     Only used for G3PLM, ignored in all other models.
#' @param M list(min, max); named list of two integer values for the minimum and maximum values for a discrete uniform draw of number of categories per item. Ignored or G3PLM.
#' @return Itembank. A list with class MCAT.items of parameters simulated under the given model.
#'     Additionally the itembank includes a list with options used to create it.
#' @export
genItembank <- function(Q=2,K=50,model='GPCM',covar=diag(Q),a=list(method="normal",a=0,b=1),b=list(method="normal",a=0,b=1),c=list(a=0,b=0),M=list(min=2,max=4)){
  # TODO: Check alpha, how to sample multivariate lognormal?
  # TODO: Check what values make sense for sampling beta, incorporate M list into beta list.
  # TODO: Check input, generate errors and/or default input.
  
  # start output list.
  out <- list()
  
  # alphas
  alpha <- matrix(rnorm(Q*K),ncol=Q)
  # correlate (eigen decomposition)
  lambda <- with(eigen(covar), vectors %*% diag(sqrt(values)))
  alpha <- t(lambda %*% t(alpha))
  # make lognormal and scale sum(alpha) to 1.
  # NOTE: Is the change in covariance of alpha important?
  alpha <- exp(alpha)
  # NOTE: What is the effect of scaling? At the very least with Q=2 it gives a perfect correlation...
  # alpha <- t(apply(alpha,1,function(x) x/sum(x)))
  
  out$alpha <- alpha
  
  #beta/eta/m
  m <- integer(K)
  if(model=="G3PLM") M <- list(min=1,max=1)
  pars <- matrix(NA,ncol=M$max,nrow=K)
  if(model=="GPCM") pars2 <- pars
  for (i in 1:K){
    m[i] <- floor(runif(1)*(M$max-M$min+1))+M$min
    par <- rnorm(m[i])
    
    # GRM requires monotously increasing beta
    if(model=="GRM" | model=="GPCM") par <- sort(par)
    pars[i,1:m[i]] <- par
    
    # for GPCM beta is actually a 'rolling sum' over eta, the pars.
    if(model=="GPCM"){
      par2 <- par
      for (j in seq_along(par)){
        par2[j] <- sum(par[1:j])
      }
      pars2[i,1:m[i]] <- par2
    }
  }
  
  # guessing (because why not?), in G3PLM only.
  if (model=="G3PLM"){
    out$guessing <- runif(K,c$a,c$b)
  }
  
  if (model=="GPCM"){
    out$eta <- pars
    out$beta <- pars2
  } else {
    out$beta <- pars
  }
  
  out$m <- m
  
  # set class and pass along options
  attr(out,"class") <- "MCAT.items"
  out$options <- list(Q=Q,K=K,M=M,model=model,covar=covar,c=c)
  out
}

#' Print itembank with some useful detail.
print.MCAT.items <- function(itembank){
  # TODO: Actually do something useful....
  for (part in itembank) print(part)
}

#` Plot itembank with some useful detail.
plot.MCAT.items <- function(itembank){
  # TODO: Actually do something useful....
  if(!require(aplpack)) install.packages('aplpack'); require(aplpack)
  plot(faces(itembank$beta))
}

subset.MCAT.items <- function(itembank,ss){
  out <- list()
  nom <- names(itembank)
  for (i in 1:length(itembank)){
    if(nom[i] != "options"){
      if(!is.vector(itembank[[i]])) out[[nom[i]]] <- matrix(itembank[[i]][ss,],ncol=ncol(itembank[[i]]))
      else out[[nom[i]]] <- itembank[[i]][ss]
    }
  }
  out$options <- itembank$options
  out$options$subset <- ss
  out
}