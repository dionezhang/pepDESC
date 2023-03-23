#' Title Normalization
#' This optinal functions performs median normalization among samples
#'
#' @param Peptide Quantitative peptide without untrusted data and contaminant data, and without outlier samples
#' @param method "Median" or "Mean", indicating to normalize by median value or mean value.
#' @param group T or F, normalizing the peptide abundances within each group or between two sample groups.
#'
#' @return data list of normalized peptide information
#' @examples
#' Peptide<-Initiation(test_raw)
#' Peptide<-FilterData(Peptide,type=test_type,'A','B')
#' Peptide<-FilterSample(Peptide)
#' Peptide<-Normalization(Peptide)
Normalization<-function(Peptide,method='Median',group=FALSE){
  s<-(min(apply(Peptide$Abu,2,median))==0)
  Norm<-function(Peptide){
    if(method=='Median'){Median<-apply(Peptide$Abu,2,stats::median)}
    if(method=='Mean'|s){Median<-apply(Peptide$Abu,2,function(x) mean(x[x>0]))}
    Peptide$Abu <- sweep(Peptide$Abu, 2, Median/mean(Median),
                         "/")
    return(Peptide$Abu)
  }
  Subset<-function(Peptide,x){
    Peptide$Abu<-Peptide$Abu[,Peptide$type==x]
    Peptide$type<-Peptide$type[Peptide$type==x]
    return(Peptide)
  }
  if(!group){
    Peptide$Abu<-Norm(Peptide)
  }
  if(group){
    for( i in unique(Peptide$type)){
      Peptide$Abu[,Peptide$type==i]<-Norm(Subset(Peptide,i))
    }
  }
  return(Peptide)
}
