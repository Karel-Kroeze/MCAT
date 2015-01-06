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
#' @export
info <- function(theta, items, method="FI", model=items$model, prior=items$prior){
  # TODO: check input.
  if(method == "FI" | method == "PFI"){
    p <- prob(theta,items,model)
    a <- items$alpha
    # Maximize determinant of Fisher Information (Segall, 1996, 2000)
    if(model == "3PL"){
      c <- items$guessing
      q <- p[,1]; p <- p[,2]
      D <- (q/p) * ((p-c)/(1-c))^2
    } else {
      K <- items$K
      m <- items$m
      D <- numeric(K)
      # straight out copy from deriv
      if (model == "GRM" | model == "SM"){
        # TODO: See if we should put the whole d1, d2 thing in prob. We're recalculating the building blocks here...
        lf <- function(x){ exp(x)/(1+exp(x)) } # from prob.
        at <- apply(a * drop(theta),1,sum)
        b <- items$beta
        
        if (model == "GRM"){
          # Graded Response Model (Glas & Dagohoy, 2007)
          for(i in 1:K){
            for(j in 1:(m[i]+1)){
              Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
              D[i] <- D[i] + p[i,j] * (Psi[j] * (1-Psi[j]) + Psi[j+1] * (1-Psi[j+1]))
            }
          }
        } else {
          # Sequential Model (Tutz, xxxx)
          # TODO: Wording in Glas & Dagohoy for dij (SM) is odd. So Psi_0 = 1, right? Why 1-Psi_(m+1)=1 and not Psi_(m+1)=0?
          # TODO: also, why is there an extra set of parenthesis in the formula for dij?
          # TODO: Finally, it doesn't work. See comments in estiamte.R.
          for(i in 1:K){
            for(j in 1:(m[i]+1)){
              Psi <- c(1,lf(at[i]-b[i,1:m[i]]),0)
              D[i] <- D[i] + p[i,j] * sum(Psi[2:(j+1)] * (1 - Psi[2:(j+1)]))
            }
          }
        }
      }
      
      if (model == "GPCM"){
        # Generalized Partial Credit Model (Muraki, 1992)
        # TODO: again, why the extra parentheses? 
        for(i in 1:K){
          mi <- 1:m[i]
          for(j in 1:(m[i]+1)){
            mp <- sum(mi*p[i,j])
            D[i] <- D[i] + p[i,j] * sum((mi * p[i,j]) * (mi - mp))
          }
        }
      }
    }
    
    out <- matrix(0,items$Q,items$Q)
    # TODO: Fisher Information for other models.
    for (i in 1:nrow(a)) { # check this for available items
      out <- out + (a[i,] %*% t(a[i,])) * D[i]
    }
    if(model == "PFI") out <- out + solve(items$prior)
  }
  
  return(out)
}

next.item <- function(theta,items,responded,available=(1:items$K)[-responded],method="FI",model=items$model,prior=items$prior,debug=FALSE) {
  I <- as.numeric(rep(NA,items$K))
  if (method == "FI" | method == "PFI"){
    # info so far
    I0 <- info(theta,subset(items,responded),method,model,prior)
    # calc det((P)FI) for available items.
    for (i in available) {
      I[i] <- det(info(theta,subset(items,i),"FI",model,prior) + I0) # FI because prior is included in I0
    }
    #print(I)
  }
  
  if(debug) print(available)
  if(debug) print(responded)
  k <- which(I == max(I,na.rm=TRUE))
  if (k < 0) k <- sample(available,1) # TODO: this is a quick hack to always return something, I'm sure we can do better.
  return(k)
}