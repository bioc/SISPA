#'@name freqplot
#'@aliases freqplot
#'@title A plotting function for SISPA sample identifiers
#'@description Given a sample changepoint data frame, will plot number of samples with and without profile activity
#'@usage freqplot(x)
#'@param x : A data frame containing samples as rows followed by zscores and estimated changepoints to be plotted.
#'@details This function expects the output from cptSamples function of SISPA package, and shows the number of samples with (orange filled bars) and without profile activity (grey filled bars).
#'@return Bar plot pdf illustrating distribution of samples
#'@import plyr
#'@examples
#'samples <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10")
#'zscores <- c(3.83,2.70,2.67,2.31,1.70,1.25,-0.42,-1.01,-2.43,-3.37)
#'changepoints <- c(1,1,1,2,2,3,3,NA,NA,NA)
#'sample_groups <- c(1,1,1,0,0,0,0,0,0,0)
#'my.data = data.frame(samples,zscores,changepoints,sample_groups)
#'freqplot(my.data)
#'@export
freqplot = function(x){
    count_data <- count(x,"sample_groups")
    freq <- barplot(count_data$freq,col=c("grey","orange"),ylab="# of Samples",names.arg=c("Inactive Profile","Active Profile"),cex.names=1.0)
    freq <- box()
    return(freq)
}


