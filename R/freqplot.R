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
#'g <- 50 ## number of genes
#'s <- 50 ## number of samples
#'## sample data matrix with values ranging from 1 to 50
#'expr <- matrix(sample.int(50, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(paste("g", 1:s, sep="") , paste("s", 1:g, sep="")))
#'## genes of interest
#'genes <- list(set1=paste("g", 1:3, sep=""))
#'## Estimates GSVA enrichment zscores.
#'gsva_results <- callGSVA(expr,genes)
#'cpt_on_samples <- cptSamples(gsva_results,dir="up",cpt_data="var",cpt_method="BinSeg",cpt_max=60)
#'## Plot number of samples by sample groups
#'freqplot(cpt_on_samples)
#'@export
freqplot = function(x){
    count_data <- count(x,"sample_groups")
    freq <- barplot(count_data$freq,col=c("grey","orange"),ylab="# of Samples",names.arg=c("Inactive Profile","Active Profile"),cex.names=1.0)
    freq <- box()
    return(freq)
}


