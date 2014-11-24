#### internal helper functions.
#' Inner product
#' 
#' Calculates the dot or inner product of any two vectors or scalars.
#' 
#' R is extremely fiddly with matrix dimensions, 
#' creating explicit functions for the desired operations will hopefully alleviate this mess slightly.
#' For inp(), order is irrelevant
#' 
#' @param x,y Vector or scalar
#' @value scalar The inner product of x and y.
inp <- function(x,y){
  x <- matrix(x,nrow=1)
  y <- matrix(y,ncol=1)
  print(x)
  print(y)
  return(drop(x %*% y))
}

#inp(1:3,2)
