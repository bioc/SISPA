#'@name sortData
#'@aliases sortData
#'@title Sorts the data by a column
#'@description Sorts the data frame by a column index in the given order 
#'@usage sortData(x,i,b)
#'@param x A data frame 
#'@param i A numeric column index of the data frame to sort it by
#'@param b User specified sorting order, ascending (FALSE) or descending (TRUE)
#'@details defaults are used: i = 1, b = FALSE, if not specified
#'@return sorted data by the input column index
#'@author Bhakti Dwivedi & Jeanne Kowalski
#'@examples
#'samples <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10")
#'zscores <- c(3.83,2.70,2.67,2.31,1.70,1.25,-0.42,-1.01,-2.43,-3.37)
#'my.data = data.frame(samples,zscores)
#'sortData(my.data,2,TRUE)
#'@export
sortData = function(x,i,b) {
  if(missing(x)){
    stop("input data is missing!")
  }
  if(missing(b)){
    b=FALSE
  }
  if(missing(i)){
    i=1
  }
  return( x[ order(x[,i], decreasing=b), ] )
}