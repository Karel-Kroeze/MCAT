#' Simulate an MCAT for a single person.
#' 
#' Wrapper to simulate an MCAT with the given criteria for a given person and itembank.
#' 
#' DETAILS (note that only PFI/FI are implemented, and only for 3PL)
#' 
#' @param person List of person parameters and related information. See \code{\link{MCAT::person}}.
#' @param itembank List of item parameters and related information. See \code{\link{MCAT::genItembank}}.
#' @param criteria List of test criteria as defined by \code{\link{MCAT::MCAT}}.
#' @return person
#' @export
MCAT <- function(person, items, st = 5, le = 50, criteria = NULL,debug=FALSE){
  if (is.null(criteria)){ # very simple defaults; TODO: create proper model check helper function. (also for items / person).
    criteria$model <- items$model
    criteria$estimate$method <- "BM"
    criteria$select$method <- "PFI"
    criteria$prior <- diag(items$Q)
  }
  process <- list()
  # very basic start.
  # set up response and available vectors.
  person$resp <- rep(NA,items$K) 
  # answer st random items to get us started.
  start <- sample(1:items$K,st)
  person$resp[start] <- ans(person$true,subset(items,start),model=criteria$model)
  person$done <- start
  
  # get theta estimate
  person$theta <- rep(0,items$Q)
  person$theta <- est(person$theta,
                      subset(items,person$done),
                      resp=person$resp[person$done],
                      method=criteria$estimate$method,
                      model=criteria$model,
                      prior=criteria$prior)
  # TODO: proper variance, this is a hacked BM variance.
  person$var <- diag(info(person$theta,subset(items,person$done),'PFI')^(-1))
  
  
  # do 25 items (total of 30)
  person$deviation <- person$theta - person$true
  if (st >= le) return(person)
  for (i in st:le){
    k <- next.item(person$theta,items,person$done,method=criteria$select$method,model=criteria$model,prior=criteria$prior,debug=debug)
    #cat("cycle:",i,"\n")
    #cat("next item:",k,"\n")
    #cat("theta:",person$theta,"\n")
    #cat("true:",person$true,"\n")
    person$resp[k] <- ans(person$true,subset(items,k),model=criteria$model)
    person$done <- c(person$done,k)
    person$theta <- est(person$theta,subset(items,person$done),resp=person$resp[person$done],method=criteria$estimate$method,model=criteria$model,prior=criteria$prior)
    # TODO: better control, but this might be able to keep ML estimates in check by limting crazyness in item selection and starting point for next cycle.
    person$theta[which(person$theta > 5)] <- 5
    person$theta[which(person$theta < -5)] <- -5
    
    # TODO: proper variance, this is a hacked BM variance.
    person$var <- diag(info(person$theta,subset(items,person$done),'PFI')^(-1))
  }  
  if(debug) print(process)
  person$deviation <- person$theta - person$true
  return(person)
}
#items <- genItembank(Q=2,model='3PL',K=100)
#MCAT(list(true=rnorm(2)),items,debug=FALSE)
