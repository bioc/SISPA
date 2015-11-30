#'@name freqplot
#'@aliases freqplot
#'@title A plotting function for SISPA sample identifiers
#'@description Given a sample changepoint data frame, will plot number of samples with and without profile activity
#'@usage freqplot(x)
#'@param x A data frame containing samples as rows followed by zscores and estimated changepoints to be plotted.
#'@details This function expects the output from cptSamples function of SISPA package, and shows the number of samples with (orange filled bars) and without profile activity (grey filled bars).
#'@return Bar plot pdf illustrating distribution of samples
#'@import plyr
#'@import ggplot2
#'@examples
#'samples <- c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10")
#'zscores <- c(3.83,2.70,2.67,2.31,1.70,1.25,-0.42,-1.01,-2.43,-3.37)
#'changepoints <- c(1,1,1,2,2,3,3,NA,NA,NA)
#'sample_groups <- c(1,1,1,0,0,0,0,0,0,0)
#'my.data = data.frame(samples,zscores,changepoints,sample_groups)
#'freqplot(my.data)
#'@export
freqplot = function(x){
    count_data <- count(x,vars="sample_groups")
    count_data$sample_groups[count_data$sample_groups=="1"] <- "Active Profile"
    count_data$sample_groups[count_data$sample_groups=="0"] <- "Inactive Profile"
    DF1 <- melt(count_data, id.var="sample_groups")
    bp <- ggplot(DF1,  aes(x=sample_groups,y=value,fill=sample_groups))
    bp <- bp + geom_bar(stat="identity",position="stack",width=1.0, colour="white")
    bp <- bp + scale_fill_manual(breaks = c("Inactive Profile","Active Profile"), values = c("orange", "grey"))
    bp <- bp + theme(axis.title.y = element_text(colour="black",size=12.0))
    bp <- bp + theme(panel.background = element_rect(fill="NA"))
    bp <- bp + theme(panel.border = element_rect(colour = "black", fill="NA"))
    bp <- bp + theme(panel.grid.major.y = element_line(colour="NA"))
    bp <- bp + theme(panel.grid.minor = element_line(colour = "NA"))
    bp <- bp + theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5, colour="black", size=12.0))
    bp <- bp + theme(axis.text.y = element_text(hjust = 1.0, vjust = 0.5, colour="black", size=10.0))
    bp <- bp + ylab("# of Samples\n")
    bp <- bp + xlab("")
    bp <- bp + theme(legend.position = "none")
    return(bp)
}


