#'@name waterfallplot
#'@aliases waterfallplot
#'@title A plotting function for SISPA sample identifiers
#'@description Given a sample changepoint data frame, will plot all samples zscores from that data.
#'@usage waterfallplot(x)
#'@param x : A data frame containing samples as rows followed by zscores and estimated changepoints to be plotted.
#'@details This function expects the output from cptSamples function of SISPA package, and highlights the sample profile of interest in the changepoint 1 with orange filled bars.
#'@return Bar plot pdf illustrating distinct SISPA sample profiles.
#'@import ggplot2
#'@examples
#'g <- 10 ## number of genes
#'s <- 30 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(paste("g", 1:g, sep="") , paste("s", 1:s, sep="")))
#'## genes of interest
#'genes <- list(set1=paste("g", 1:3, sep=""))
#'## Estimates GSVA enrichment zscores.
#'gsva_results <- callGSVA(expr,genes)
#'cpt_on_samples <- cptSamples(gsva_results,dir="up",cpt_data="var",cpt_method="BinSeg",cpt_max=60)
#'waterfallplot(cpt_on_samples)
#'@export
waterfallplot = function(x){
    if(missing(x)){
      stop("input data is missing!")
  }
    x <- x[,-ncol(x)]
    colnames(x)[1] <- "samples"
    #create a start column
    x[,ncol(x)+1] <- 0
    #define global variable
    #select and modify the data for readability
    read_data_sort <- x[ order(x[,2], decreasing=FALSE), ]
    read_data_sort$changepoints[is.na(read_data_sort$changepoints)] <- 1000
    #arrange the columns
    arrange_read_data_sort <- data.frame(read_data_sort[,1],read_data_sort[,4],read_data_sort[,2],read_data_sort[,3])
    colnames(arrange_read_data_sort)<-c("samples","Start","End","changepoints")
    arrange_read_data_sort[,1] <- factor(arrange_read_data_sort[,1], levels=arrange_read_data_sort[,1])
    arrange_read_data_sort$id <- seq_along(arrange_read_data_sort$End)
    arrange_read_data_sort$type <- ifelse(arrange_read_data_sort$changepoints==1, "cpt1", "other")
    arrange_read_data_sort <- arrange_read_data_sort[,c(5,1,6,2,3)]
    #generate the waterfall plot
    wf <- ggplot(arrange_read_data_sort, aes(colour=type,samples,fill=type)) + scale_fill_manual(values=c("#FF8000","white")) + scale_colour_manual(values=c("#FF8000","grey")) + geom_rect(aes(x=samples,xmin=id-0.35,xmax=id+0.35,ymin=End,ymax=Start))  
    wf <- wf + labs(x="", y="", fill="")
    wf <- wf + theme(legend.position = "none")
    wf <- wf + theme(axis.text.y = element_text(colour="black",size=12.0))
    wf <- wf + theme(axis.text.x = element_blank())
    #wf <- wf + theme(axis.title.y = element_text(colour="black",face="bold",size=12.0))
    wf <- wf + theme(panel.background = element_rect(fill="NA"))
    wf <- wf + theme(panel.border = element_rect(colour = "black", fill="NA"))
    wf <- wf + theme(panel.grid.major.y = element_line(colour="NA"))
    wf <- wf + theme(panel.grid.minor = element_line(colour = "NA"))
    wf <- wf + theme(axis.ticks.x = element_blank())
    wf <- wf + theme(axis.ticks.y=element_line(colour="black"))
    wf <- wf + ylab("SISPA Profile Score")
    wf <- wf + xlab("Samples\n")
    return(wf)
}



