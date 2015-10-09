#'@name callGSVA
#'@aliases callGSVA
#'@title GSVA enrichment analysis
#'@description Estimates GSVA enrichment zscores.
#'@usage callGSVA(x,y)
#'@param x : A data matrix of gene or probe expression values where rows corrospond to genes and columns corrospond to samples
#'@param y : Gene sets provided as a list object.
#'@details This function uses "zscore" gene-set enrichment method in the estimation of gene-set enrichment scores per sample.
#'@return A gene-set by sample matrix of GSVA enrichment zscores.
#'@import GSVA
#'@examples 
#'g <- 10 ## number of genes
#'s <- 30 ## number of samples
#'## sample data matrix with values ranging from 1 to 10
#'expr <- matrix(sample.int(10, size = g*s, replace = TRUE), nrow=g, ncol=s, dimnames=list(paste("g", 1:g, sep="") , paste("s", 1:s, sep="")))
#'## genes of interest
#'genes <- list(set1=paste("g", 1:3, sep=""))
#'## Estimates GSVA enrichment zscores.
#'callGSVA(expr,genes)
#'@seealso GSVA
#'@export
callGSVA = function(x,y) {
    if(missing(x)){
    stop("input expression data missing!")
    }
    if(missing(y)){
    stop("input gene set missing!")
    }
    #rnaseq=TRUE does not work with method="zscore"
    gsva.results <- gsva(x, y, method="zscore", rnaseq=FALSE,verbose=FALSE,parallel.sz=2)
    tr_gsva.results <- t(gsva.results)
    #label column names
    tr_result_zscore <- cbind(samples = rownames(tr_gsva.results), tr_gsva.results)
    colnames(tr_result_zscore) <- c("samples","zscore")
    rownames(tr_result_zscore) <- NULL
    return ( as.data.frame(tr_result_zscore))
}