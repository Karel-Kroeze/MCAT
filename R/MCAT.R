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
MCAT <- function(person, items, criteria = NULL,debug=FALSE){
  if (is.null(criteria)){ # very simple defaults; TODO: create proper model check helper function. (also for items / person).
    criteria$model <- "3PL"
    criteria$estimate$method <- "EAP"
    criteria$select$method <- "PFI"
    criteria$prior <- diag(items$Q)
  }
  process <- list()
  # very basic start.
  # set up response and available vectors.
  person$resp <- rep(NA,items$K) 
  # answer 5 random items to get us started.
  start <- sample(1:items$K,5)
  person$resp[start] <- ans(person$true,subset(items,start),model=criteria$model)
  person$done <- start
  
  # get theta estimate
  person$theta <- rep(0,items$Q)
  person$theta <- est(person$theta,subset(items,person$done),resp=person$resp[person$done],method=criteria$estimate$method,model=criteria$model,prior=criteria$prior)

  # do 25 items (total of 30)
  for (i in 1:50){
    k <- next.item(person$theta,items,person$done,method=criteria$select$method,model=criteria$model,prior=criteria$prior,debug=debug)
    person$resp[k] <- ans(person$true,subset(items,k),model=criteria$model)
    person$done <- c(person$done,k)
    person$theta <- est(person$theta,subset(items,person$done),resp=person$resp[person$done],method=criteria$estimate$method,model=criteria$model,prior=criteria$prior)
    
    process[[i]] <- list(item = subset(items,k), resp = person$resp[k], theta = person$theta)
  }  
  if(debug) print(process)
  return(person)
}

#items <- genItembank(Q=2,model='3PL',K=100)
#MCAT(list(true=rnorm(2)),items,debug=FALSE)
