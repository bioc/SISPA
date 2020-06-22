#'@name filterVars
#'@aliases filterVars
#'@title A filter function for the data
#'@description Filter rows with zero values
#'@usage filterVars(x,y)
#'@param x : A data frame or matrix where rows represent gene and columns represent samples
#'@param y : A vector of a sample column values to apply the filtering on.
#'@details This function filter out rows with zero data value for a given sample. Both input arguments (x and y) must be of the same length
#'@return The returned value is a list containing an entry for each row filtered out by zero data value
#'@examples 
#'x = matrix(runif(3*10, 0, 1), ncol=3)
#'y <- x[,1]
#'filterVars(x,y)
#'@export
filterVars = function(x,y) {
  if(missing(x)){
    stop("input data is missing!")
  }
  if(missing(y) | !is.vector(y)){
    stop("data needs to be a vector!")
  }
  return( x[ y>"0", ] )
}