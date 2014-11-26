############# Information functions and next.item wrapper
#' Item information
#' 
#' Computes item information for a given item and theta.
#' 
#' DETAILS
#' 
#' @param theta Numeric (vector) of latent abilitities.
#' @param items Itembank
#' @param model IRT model to be used.
#' @param method Information function to be used.
#' @return information Numeric (vector) of information functions evaluated for the given item(s) and theta.
info <- function(theta, items, method="FI", model=items$model, prior=items$prior){
  # TODO: check input.
  if(method == "FI" | method == "PFI"){
    # Maximize determinant of Fisher Information (Segall, 1996, 2000)
    if(model == "3PL"){
      FI <- function(theta,items){
        out <- matrix(0,items$Q,items$Q)
        a <- items$alpha
        c <- items$guessing
        p <- prob(theta,items,model="3PL"); q <- p[,1]; p <- p[,2]
        w <- (q/p) * ((p-c)/(1-c))^2
        for (i in 1:nrow(a)) { # check this for available items
          out <- out + (a[i,] %*% t(a[i,])) * w[i]
        }
        return(out)
      }
    }
    # TODO: Fisher Information for other models.
    out <- FI(theta,items)
    if(model == "PFI") out <- out + solve(items$prior)
  }
  
  return(out)
}

next.item <- function(theta,items,responded,available=(1:items$K)[-responded],method="FI",model=items$model,prior=items$prior,debug=FALSE) {
  I <- numeric(items$K)
  
  if (method == "FI" | method == "PFI"){
    # info so far
    I0 <- info(theta,subset(items,responded),method,model,prior)
    # calc det((P)FI) for available items.
    for (i in available) {
      I[i] <- det(info(theta,subset(items,i),"FI",model,prior) + I0) # FI because prior is included in I0
    }
  }
  
  if(debug) print(available)
  if(debug) print(responded)
  return(which(I == max(I)))
}