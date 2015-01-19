#' Generate itembank.
#' 
#' Generates a Q-dimensional itembank under the Graded Response Model (GRM), Sequential Model (SM), Generalized Partial Credit Model (GPCM) or Generalized Three Parameter Logistic Model (3PL).
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
#'     Only used for 3PL, ignored in all other models.
#' @param M list(min, max); named list of two integer values for the minimum and maximum values for a discrete uniform draw of number of categories per item. Forced to 1,1 in 3PL.
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
  #alpha <- matrix(rlnorm(Q*K,0,.25),ncol=Q)
  require(mvtnorm)
  alpha <- rmvnorm(K,rep(1,Q),covar)
  # TODO; A way to correlate alpha without making it negative...
  # correlate (eigen decomposition)
  #lambda <- with(eigen(covar), vectors %*% diag(sqrt(values)))
  #alpha <- t(lambda %*% t(alpha))
  # make lognormal and scale sum(alpha) to 1.
  # NOTE: Is the change in covariance of alpha important?
  # NOTE: What is the effect of scaling? At the very least with Q=2 it gives a perfect correlation...
  # alpha <- t(apply(alpha,1,function(x) x/sum(x)))
  alpha[which(alpha < 0)] <- 0
  out$alpha <- alpha
  
  #beta/eta/m
  m <- integer(K)
  if(model=="3PL") M <- list(min=1,max=1)
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
  if (model=="3PL"){
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
  options <- list(Q=Q,K=K,M=M,model=model,covar=covar,c=c)
  out <- c(out,options)
  attr(out,"class") <- "MCAT_itembank"
  return(invisible(out))
}

#' Print itembank with some useful detail.
#' @export
print.MCAT_itembank <- function(itembank){
  # TODO: Actually do something useful....
  for (i in seq_along(itembank)){
    cat(names(itembank)[i],"\n")
    print(itembank[[i]])
  } 
}

#' Plot itembank with some useful detail.
#' @export
plot.MCAT_itembank <- function(itembank){
  # TODO: Actually do something useful....
  if(!require(aplpack)) install.packages('aplpack'); require(aplpack)
  plot(faces(itembank$beta))
}

#' Create a subset of the itembank containing the given items.
#' @param itembank MCAT_items, list with class MCAT_items, see genItembank
#' @param ss, vector of integers, indeces of item(s) to be included in the subset.
#' @return MCAT_itembank, subset containing item(s) ss and all the other bits and pieces that make an itembank.
#' @export
subset.MCAT_itembank <- function(itembank,ss){
  # TODO: fix and properly handle double subsets.
  out <- list()
  nom <- names(itembank)
  pars <- c("alpha","beta",'eta','guessing','m') # item parameters, everything else is fluff. 
  for (i in 1:length(itembank)){
    if(nom[i] %in% pars){
      if(!is.vector(itembank[[i]])) out[[nom[i]]] <- matrix(itembank[[i]][ss,],ncol=ncol(itembank[[i]]))
      else out[[nom[i]]] <- itembank[[i]][ss]
    } else {
      out[[nom[i]]] <- itembank[[i]]
    }
  }
  out$subset <- ss
  attr(out,"class") <- "MCAT_itembank"
  return(out)
}

## TODO: think on how to deal with nested subset calls, currently this fails massively (indices incorrect)