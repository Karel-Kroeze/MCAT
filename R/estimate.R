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
#' @param model String, valid values are "3PL", "GRM", "SM" or "GPCM". Defaults to "GPCM".
#' @param Sigma0 Numeric (matrix), prior covariance. Only required when method != "ML".
#' @param debug Logical. Should debug information be shown?
#' @param ... Passed on to eapEst, ignored if method is not "EAP". See \code{\link{eapEst}} for details.
#' @return list with two entries; theta contains the current theta estimate, and SE the standard error of the estimate.
#' @export
est <- function(theta, items, resp, 
                method = "BM", model = items$model,
                prior = items$covar, threshold = 1e-4, max.iter = 10,
                debug = FALSE, ip = 10,...){
  if (method == "ML") prior <- NULL
  if (method == "MAP") method <- "BM"
  
  if (method == "EAP") out <- eapEst(items,resp,model,prior,ip,debug,...)
  else {
    iter <- 0
    delta <- 1
    # Newton-Rhapson approximation.
    out <- tryCatch({
      theta.last <- theta
      while(mean(abs(delta)) > threshold & iter <= max.iter){
        theta0 <- theta
        theta <- with(deriv(theta,items,resp,model,prior),
                      theta0 - solve(d2) %*% d1)
        delta <- theta - theta0
        iter <- iter + 1
        if(debug) cat('N-R>',iter,':',theta,"(",delta,")\n")
      }
      theta
    }, warning = function(w) {
      cat(geterrmessage())
      return(theta.last)
    }, error = function(e) {
      cat(geterrmessage())
      return(theta.last)
    })
  }
  return(as.vector(out))
}