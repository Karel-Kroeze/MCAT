#' Obtain an ability estimate given a response pattern and item parameters. 
#' 
#' Calculates either the maximum likelihood (ML), estimated a-posteriori (EAP) or Bayes modal (BM or MAP) estimate.
#' 
#' Details
#' 
#' @param theta Numeric (vector), initial theta estimate. Defaults to a null vector of length equal to the number of dimensions.
#' @param items Itembank, see \code{\link{genItembank}}.
#' @param resp Response pattern, see \code{\link{ans}}.
#' @param method String, valid values are "ML", "EAP", and "MAP" or "BM". Defaults to "BM".
#' @param model String, valid values are "G3PLM", "GRM", "SM" or "GPCM". Defaults to "GPCM".
#' @param Sigma0 Numeric (matrix), prior covariance. Only required when method != "ML".
#' @param debug Logical. Should debug information be shown?
#' @param ... Passed on to eapEst, ignored if method is not "EAP". See \code{\link{eapEst}} for details.
#' @return list with two entries; theta contains the current theta estimate, and SE the standard error of the estimate.
#' @export
est <- function(theta, items, resp, 
                method = "BM", model = "GPCM",
                prior = NULL, threshold = 1e-4, max.iter = 10,
                debug = FALSE, ip = 10,...){
  if (method == "ML") prior <- NULL
  if (method == "MAP") method <- "BM"
  
  if (method == "EAP") out <- eapEst(items,resp,model,prior,ip,debug,...)
  else {
    iter <- 0
    delta <- 1
    # Newton-Rhapson approximation.
    while(mean(abs(delta)) > threshold & iter <= max.iter){
      theta0 <- theta
      theta <- with(deriv(theta,items,resp,model,prior),
                    theta0 - solve(d2) %*% d1)
      delta <- theta - theta0
      iter <- iter + 1
      if(debug) cat('N-R>',iter,':',theta,"(",delta,")\n")
    }
    out <- theta
  }
  return(as.vector(out))
}

############# TESTING CODE>
# SM "works" but only for theta <= 0, and underestimates theta. Positive theta explodes. So really, it doesn't work.
# GPCM is close to a decimal point, but differs slightly.
# GRM and G3PLM (4PL) are spot-on.
# Both ML and BM (do not) work in the same situations. 
#
# EAP works in theory, but has some significant numerical issues. At higher numbers of items the estimate approaches 1.157, not sure why that number.
# Annoyingly, it is never the same as catR. catR's numerical integration is very simple, but it isn't THAT far off...
#### Test estimators.
# Q <- 1
# theta <- rep(0,Q)
# model <- "GRM"
# method <- "BM"
# fix.alpha <- T # Set all alpha to 1?
# show.newton <- T # Show Newton-Rhapson iteration steps?
# K <- 25
# 
# max.iter <- 1000
# threshold <- 1e-4
# xmodel <- switch(model,"G3PLM"=NULL,model)
# items <- genItembank(model=model,K=K,Q=Q,)
# if(fix.alpha) items$alpha <- matrix(1,K,Q)
# resp <- ans(theta,items,model)
# prior <- diag(Q) 
# xitems <- switch(model,
#                  "G3PLM"=cbind(items$alpha,items$beta,items$guessing,1), # 4PLM with upper asymptote = 1
#                  "GPCM"=cbind(items$alpha,items$eta),
#                  cbind(items$alpha,items$beta))
# 
# cat(method," estimates for ",model," model (",K," items):\n",sep="")
# cat("TRUE:",theta,"\n")
# cat("MCAT:",est(theta,items,resp,method,model,prior,threshold,max.iter,debug=show.newton),"\n")
# if(Q == 1) cat("catR:",thetaEst(xitems,resp,xmodel,1,method),"\n")
# 

### timings. Ouch, catR is up to 10 times faster...
# TODO: Look into a native newton-rhapson or something similar. Rootfinding sadly won't work?
# system.time(replicate(10,est(theta,items,resp,method,model,prior,threshold,max.iter)))
# system.time(replicate(10,thetaEst(xitems,resp,xmodel,1,method)))
