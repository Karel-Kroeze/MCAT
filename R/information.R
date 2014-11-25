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
#' @value information Numeric (vector) of information functions evaluated for the given item(s) and theta.
info <- function(theta, items, model, method){
  # TODO: check input.
  if(method == "FI"){
    # Maximize determinant of Fisher Information (Segall, 1996, 2000)
    if(model == "3PL"){
      
    }
  }
}

info(0,0,0,"FI")
