######## MCAT test functions ########

#' Select the next item, given a test and person
#' 
#' Returns the next item index, based on item selection criteria, constraints and other applicable criteria.
#' 
#' @param MCAT test object, see \code{\link{initMCAT}}.
#' @param person MCAT_person object, see \code{\link{initPerson}}
#' @param debug boolean, print debug info?
#' @export
next_item <- function(MCAT,person,debug=FALSE) {
  if (MCAT$start$type == "random"){
    if (MCAT$start$n > length(person$done)){
      return(sample(person$avail,1))
    }
  } else if (MCAT$start$type == "fixed"){
    if (length(MCAT$start$items) > length(person$done)){
      print(MCAT$start)
      return(MCAT$start$items[length(person$resp) + 1])
    }
  }
  
  I <- as.numeric(rep(NA,MCAT$items$K))
  if (MCAT$itemSelection == "FI" | MCAT$itemSelection  == "PFI"){
    # info so far
    if (length(person$done) > 0){
      I0 <- info(person$theta,subset(MCAT$items,person$done),MCAT$itemSelection,MCAT$model,person$prior)
    } else {
      I0 <- matrix(0,ncol=MCAT$items$Q,nrow=MCAT$items$Q)
    }
    # calc det((P)FI) for available items.
    for (i in person$avail) {
      I[i] <- det(info(person$estimate,subset(MCAT$items,i),"FI",MCAT$model,perosn$prior) + I0) # FI because prior is included in I0
    }
  }
  
  k <- which(I == max(I,na.rm=TRUE))
  if (k < 0) k <- sample(person$avail,1) # TODO: this is a quick hack to always return something, I'm sure we can do better.
  return(k)
}

#' Update person object after a response. Update administrative variables, and calculate new estimates and variance.
#' @param MCAT MCAT test object, see \code{\link{initMCAT}}.
#' @param person MCAT_person object, see \code{\link{initPerson}}
#' @param item integer, index of item responded to.
#' @param answer integer, actual response.
#' @param ... parameters passed onto the estimator. threshold and max.iter influence BM and ML Newton-Rhapson estimation routine, ip sets the number of quadrature points per dimension in EAP estimation.
#' @return MCAT_person updated person object.
#' @export
update_response <- function(MCAT,person,item,answer,...) {
  person$avail <- person$avail[person$avail != item]
  person$done <- c(person$done,item)
  person$resp[item] <- answer
  person$estimate <- est(person$estimate,MCAT$items,person$resp,MCAT$est,MCAT$model,person$prior,...)
  person$var <- info(person$estimate,subset(MCAT$items,person$done),'PFI')^(-1) 
  # TODO: implement proper variance function, this is a hacky BM var for all types. 
  return(person)
}