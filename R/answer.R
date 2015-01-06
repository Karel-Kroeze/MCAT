#' Simulate a response (pattern) with a given theta and item(s)
#' 
#' @param theta Numeric (vector). Vector of length Q, where each element is the theta for the corresponding dimension
#' @param items Itembank, or a subset of an itembank. See \code{\link{genItembank}} and \code{\link{subset.MCAT.items}}
#' @return integer (vector) of simulated response(s).
#' @export
ans <- function(theta,items,model=items$model){
  Pij <- prob(theta = theta,items = items,model = model)
  cp <- Pij
  for (i in 1:ncol(cp)) cp[,i] <- apply(matrix(Pij[,1:i],ncol=i),1,sum)
  rand <- runif(nrow(items$beta))
  out <- apply(rand > cp,1,sum,na.rm=TRUE)
  # cat(cp,"\n",rand," -> ",out,"\n\n") # debug
  return(out)
}

# TODO: doublecheck this, it was too easy.