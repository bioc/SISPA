#'@name SISPA
#'@aliases SISPA
#'@title SISPA
#'@description SISPA: Method for Sample Integrated Gene Set Analysis
#'@usage SISPA(x,y)
#'@param x : A data matrix of gene or probe expression values where rows corrospond to genes and columns corrospond to samples
#'@param y : Gene sets provided as a list object.
#'@details Sample Integrated Gene Set Analysis (SISPA) is a method designed to define sample groups with similar gene set enrichment profiles. The user specifies a gene list of interest and sample by gene molecular data (expression, methylation, variant, or copy change data) to obtain gene set enrichment scores by each sample. The score statistics is rank ordered by the desired profile (e.g., upregulated or downregulated) for samples. A change point model is then applied to the sample scores to identify groups of samples that show similar gene set profile patterns. Samples are ranked by desired profile activity score and grouped by samples with and without profile activity. Figure 1 shows the schematic representation of the SISPA method overview.
#'@return The input molecular data frame with added sample identifiers and estimated changepoints. A plot showing the changepoint locations estimated on the data. Bar plots pdf illustrating distinct distribution of samples with and without profile activity
#'@import GSVA
#'@import changepoint
#'@import data.table
#'@import plyr
#'@import ggplot2
#'@include callGSVA.R
#'@include cptSamples.R
#'@include waterfallplot.R
#'@include freqplot.R
#'@examples
#'g <- 10 ## number of genes
#'s <- 30 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(paste("g", 1:g, sep="") , paste("s", 1:s, sep="")))
#'## genes of interest
#'genes <- list(set1=paste("g", 1:3, sep=""))
#'SISPA(expr,genes)
#'@export
SISPA = function(x,y){
  
#sample enrichment scores: GSVA
    gsva_results <- callGSVA(x,y)

#sample profile identification
    ##pdf(file="changepoints_all.pdf",height=6,width=6)
    cpt_on_samples <- cptSamples(gsva_results,dir="up",cpt_data="var",cpt_method="BinSeg",cpt_max=60)
    ##dev.off()
    write.table(cpt_on_samples,file="changepoints_all.csv",row.names=FALSE,sep=",")

#sample profile distribution
    #waterfall plot
    ##pdf(file="waterfallplot.pdf",height=6,width=6)
    wfplot <- waterfallplot(cpt_on_samples)
    ##print(wfplot)
    ##dev.off()
    #barplot
    ##pdf(file="freqplot.pdf",height=6,width=6)
    fqplot <- freqplot(cpt_on_samples)
    ##print()
    ##dev.off()

}