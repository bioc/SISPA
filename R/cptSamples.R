#'@name cptSamples
#'@aliases cptSamples
#'@title Sample profile identifier analysis
#'@description Generate sample profile identifiers from sample zscores using change point model.
#'@param x A matrix or data frame of sample GSVA enrichment zscores within which you wish to find a changepoint.
#'@param cpt_data Identify changepoints for data using variance (cpt.var) or mean (cpt.mean). Default is cpt.var.
#'@param cpt_method Choice of single or multiple changepoint model. Default is "BinSeg".
#'@param cpt_max The maximum number of changepoints  to search for using "BinSeg" method. Default is 60.
#'@details This function assigns samples identified in the first changepoint with the active profile ("1") while the remaining samples are grouped under inactive profile ("0").
#'@return The input data frame with added sample identifiers and estimated changepoints. A plot showing the changepoint locations estimated on the data
#'@seealso changepoint
#'@import plyr
#'@import data.table
#'@include sortData.R
#'@include callZSCORE.R
#'@examples
#'g <- 10 ## number of genes
#'s <- 60 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'rnames <- paste("g", 1:g, sep="")
#'cnames <- paste("s", 1:s, sep="")
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(rnames, cnames))
#'## genes of interest
#'genes <- data.frame(paste("g", 1:6, sep=""))
#'## Estimates GSVA enrichment zscores.
#'gsva_results <- callGSVA(expr,genes)
#'cptSamples(gsva_results,cpt_data="var",cpt_method="BinSeg",cpt_max=60)
#'@export
cptSamples <- function (x,cpt_data,cpt_method,cpt_max){
  ## Identify the changepoints for the data
  # x : data within which you wish to identify the changepoints
  # i : numeric column index of the data frame to sort it by
  # b : sorting order, ascending (FALSE) or descending (TRUE)
  perDiffcpt = function(x,cpt_data,cpt_method,cpt_max){
    max <- length(x)
    if (cpt_method == "AMOC") {
      if(cpt_data == "mean"){
        changepoints <- cpt.mean(x,method=cpt_method,Q=1)
      }
      else if(cpt_data == "var"){
        changepoints <- cpt.var(x,method=cpt_method,Q=1)
      }
    }
    else if (cpt_method == "PELT" || cpt_method == "SegNeigh" || cpt_method == "BinSeg") {
      if(max > cpt_max){
        if(cpt_data == "mean"){
          changepoints <- cpt.mean(x,method=cpt_method,Q=cpt_max)
        } else{
          changepoints <- cpt.var(x,method=cpt_method,Q=cpt_max)
        }
      }
      if(max <= cpt_max){
        if(cpt_data == "mean"){
          changepoints <- cpt.mean(x,method=cpt_method,Q=5)
        } else{
          changepoints <- cpt.var(x,method=cpt_method,Q=5)
        }
      }
      
    }
    return ( cpts(changepoints) )
  }
  ## Add the estimated changepoint locations to the data frame
  # x : A data frame containing the data used for change point estimation
  # y : vector containing the changepoint locations for the data supplied
  cptAdd = function(x,y){
    level=character()
    l = 1
    for(i in 1:length(y) ) {
      if(i==1){
        for(j in 1:y[i]){
          level[j]=l
        }
      }else{
        a <- (y[i-1])+1
        for(j in a:y[i]){
          level[j]=l
        }
      }
      if(i==length(y)){
        l=1000
        a <- y[i]+1
        for(j in a:nrow(x)){
          level[j]=1000
        }
      }
      
      l=l+1
    }
    return(as.numeric(level))
  }
  ## Locate all values for the estimated change point locations
  # x : A data frame containing the data values used for change point estimation
  # y : Vector containing the changepoint locations for the data supplied
  # z : Column index of the data values used for changepoint estimation
  cptLocate_all = function(x,y,z){
    cutoff_all <- character()  
    for(i in 1:length(y) ) {
      cutoff_all[i] <- x[y[i],][z]
    }
    cpts_plot_cutoff_all <- cutoff_all[!is.na(cutoff_all)]
    cpts_plot_cutoff_all <- as.numeric(cpts_plot_cutoff_all)
    cpts_plot_cutoff_all <- round(cpts_plot_cutoff_all,digits=2)
    return(cpts_plot_cutoff_all)
  }
  ## Plot the single or multiple changepoints identified within the data
  # x : A vector containing the data within which changepoints were identified
  # y : A vector containing the identified number of changepoints on the data
  # LABEL : Y-axis label
  cptsPlot <- function(x,y,LABEL){
    if(missing(x)){
      stop("data is missing!")
    }
    if(missing(y)){
      stop("changepoint cutoffs are missing!")
    }
    x <- sort(x,FALSE)
    y <- sort(y,TRUE)
    
    cptplot <- plot(x, type='p', lty=3, lwd=1, main="", xlab="Samples", ylab=LABEL, cex=1.5, pch=1, xlim=c(0,length(x)))
    
    for(i in 1:length(y)) {
      if(i==1) {
        abline(h=y[i],lty=2,col='red')
        text(y[i],y[i],y[i], cex=1.0, pos=4, col="red")
      }else{
        ## display in terms of sample group profiles
        if(y[i]<0){
          abline(h=y[i],lty=2,col='blue')
          text(y[i],y[i],y[i], cex=1.0, pos=4, col="blue")
        } else{
          abline(h=y[i],lty=2,col='green3')
          text(y[i],y[i],y[i], cex=1.0, pos=4, col="green3")
        }
      }
    }
    
    return(cptplot)
  }
    ######### END: sub-functions used by the cptSamples main function #########

    #read sample and sample zscores
    st <- x[,c(1,2)]
    colnames(st)[2] <- "zscore"
    #correct for empty cells or NAs
    st_rmv_na <- st[!is.na(st[ncol(st)]), ]
    #check if the zscores are numeric and sort them
    st_rmv_na$zscore <- as.numeric(as.character(st_rmv_na$zscore))
    st_rmv_na_sort = sortData(st_rmv_na,ncol(st_rmv_na),'TRUE')
    
    #identify the changepoints
    cpts = perDiffcpt(st_rmv_na_sort$zscore,cpt_data,cpt_method,cpt_max)
    #check for the number of changepoints identified
    if(length(cpts)<1){
      warning("No changepoints identified in the data set!")
      cpt_out=NULL
    }else{
      #add estimated changepoint locations to the data
      changepoints <- cptAdd(st_rmv_na_sort,cpts)
   
      #append changepoints to the original data frame
      cpt_out <- data.frame(st_rmv_na_sort,changepoints)
      all_cutoffs <- cptLocate_all(cpt_out,cpts,2)
      cpt_out_sort <- sortData(cpt_out,ncol(cpt_out)-1,'FALSE' )
  
      #plot all changepoints identified
      Y_LABEL <- "Zscores"
      ##pdf(file="changepoints_all.pdf",width=6,height=6)
      cptsPlot(cpt_out_sort$zscore,all_cutoffs,Y_LABEL)
      ##dev.off()

      #generate the data frame with all changepoints
      #format the data file for reading usability
      #any changepoint other than 1 is designated as '0'
      #plot changepoints within the data identified
      cpt_out$sample_groups <- cpt_out$changepoints
      cpt_out [, ncol(cpt_out) ][ cpt_out [ ,ncol(cpt_out) ] >1 ] <- "0"
      cpt_out [, ncol(cpt_out)-1 ][ cpt_out [ ,ncol(cpt_out)-1 ] == 1000 ] <- "NA"
    }
    return(cpt_out)
}


  