#'@name SISPA
#'@aliases SISPA
#'@title SISPA
#'@description SISPA: Method for Sample Integrated Gene Set Analysis
#'@usage SISPA(feature=1,f1.df,f1.genes,f1.profile,f2.df,f2.genes,f2.profile,cpt_data="var",cpt_method="BinSeg",cpt_max=60)
#'@param feature Number of input feature or data types
#'@param f1.df A data matrix of first feature (e.g., gene or probe expression values) where rows corrospond to genes and columns corrospond to samples
#'@param f1.genes Gene sets for first feature provided as a vector or data frame
#'@param f1.profile A flag to specify gene profile. If gene.profile="up" then samples with increased zscores are identified. If gene.profile="down" then samples with decreased zscores are identified. Default is "up".
#'@param f2.df A data matrix of second feature (e.g., gene variant change) where rows corrospond to genes and columns corrospond to samples
#'@param f2.genes Gene sets for second feature provided as a vector or data frame
#'@param f2.profile A flag to specify gene profile. If gene.profile="up" then samples with increased zscores are identified. If gene.profile="down" then samples with decreased zscores are identified. Default is "up".
#'@param cpt_data Identify changepoints for data using variance (cpt.var) or mean (cpt.mean). Default is cpt.var.
#'@param cpt_method Choice of single or multiple changepoint model. Default is "BinSeg". See changepoint R package for details
#'@param cpt_max The maximum number of changepoints  to search for using "BinSeg" method. Default is 60.
#'@details Sample Integrated Gene Set Analysis (SISPA) is a method designed to define sample groups with similar gene set enrichment profiles. The user specifies a gene list of interest and sample by gene molecular data (expression, methylation, variant, or copy change data) to obtain gene set enrichment scores by each sample. The score statistics is rank ordered by the desired profile (e.g., upregulated or downregulated) for samples. A change point model is then applied to the sample scores to identify groups of samples that show similar gene set profile patterns. Samples are ranked by desired profile activity score and grouped by samples with and without profile activity. Figure 1 shows the schematic representation of the SISPA method overview.
#'@return The input molecular data frame with added sample identifiers and estimated changepoints. A plot showing the changepoint locations estimated on the data. Bar plots pdf illustrating distinct distribution of samples with and without profile activity
#'@import GSVA
#'@import changepoint
#'@import data.table
#'@import plyr
#'@import ggplot2
#'@include sortData.R
#'@include callGSVA.R
#'@include callZSCORE.R
#'@include cptSamples.R
#'@include waterfallplot.R
#'@include freqplot.R
#'@examples
#'g <- 10 ## number of genes
#'s <- 60 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'rnames <- paste("g", 1:g, sep="")
#'cnames <- paste("s", 1:s, sep="")
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(rnames, cnames))
#'## genes of interest
#'genes <- data.frame(paste("g", 1:6, sep=""))
#'SISPA(feature=1,f1.df=expr,f1.genes=genes,f1.profile="up")
#'@export
SISPA = function(feature=1,f1.df,f1.genes,f1.profile,f2.df=NULL,f2.genes=NULL,f2.profile=NULL,cpt_data="var",cpt_method="BinSeg",cpt_max=60){
  
#when only one feature sample profile is estimated
  if(feature==1){
    if(ncol(f1.df)<4){
      stop("Zscore cannot be computed for samples < 3!")
    }
  } else {
    #zscores can not be computed if number of samples <3
    if(ncol(f1.df)<4 | ncol(f2.df)<4){
      stop("Zscore cannot be computed for samples < 3!")
    }
  }
  
#loop over each feature
for(f in 1:feature){
#sample enrichment scores: GSVA/ZSCORES
    #when number of genes is less than 3
    if (nrow(get(paste("f",f,".genes",sep="")))<3){
      #compute zscores for the gene set
      scores <- callZSCORE(get(paste("f",f,".df",sep="")))
    }else{#when number of genes more than or equal to 3
      #compute zscores using GSVA
      scores <- callGSVA(get(paste("f",f,".df",sep="")),get(paste("f",f,".genes",sep="")))
    }
#determine direction based on gene.profile
    if(get(paste("f",f,".profile",sep="")) == "up"){
      #increased zscore values are desirable
    }else{
      #decreased zscore values are desirable
      #Thus multiplied by -1 to obtain the desired outcome
      scores$NegOneValue_zscore <- scores$zscore * -1
      scores <- scores[,-c(2)]
      colnames(scores)[2] <- "zscore"
    }
    #save the feature profile scores
    assign(paste("feature",f,sep=""),scores)
    
    if(f==2){
      #sum the scores over the two features
      summed_scores <- merge(feature1,feature2,by="samples",sort=FALSE)
      summed_scores$zscores <-  as.numeric(as.character(summed_scores$zscore.x)) + as.numeric(as.character(summed_scores$zscore.y))
      summed_scores <- summed_scores[,c(1,ncol(summed_scores))]
    }else{
      summed_scores <- scores
    }
}

#sample profile identification
    cpt_on_samples <- cptSamples(summed_scores,cpt_data,cpt_method,cpt_max)

#sample profile distribution
    if(is.null(cpt_on_samples)){
      warning("No changepoints identified in the data set!")
      return (0)
    }else{
      return(cpt_on_samples)
    }
}