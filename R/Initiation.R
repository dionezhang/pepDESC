#' Title Initialization
#' This function initialize the input data to enable downstream analysis
#' @param data
#' a .Rdata file read from the peptide-level searching result
#' should contain columns "Accession" "Description" "Sequence" "Modifications"
#' Columns for peptide abundance should be annotated with "Abundance"
#' Better contain columns "RT" and "m/z"
#'
#' @return A data list of peptide information for downstream analysis
#'
#' @examples
#' Peptide<-Initiation(test_raw)
Initiation<-function (data) {
  rownames(data) <- seq_along(rownames(data))
  Check<-function(data,name){
    return(!(FALSE%in%vapply(name,function(x) x%in%colnames(data),c(1))))
  }
  stopifnot('Not enough information for initialization'=
              Check(data=data,name=c('Sequence','Modifications','Accession','Description')))
  if(!Check(data=data,name=c('RT','m/z'))){PeakInfo=FALSE}else{PeakInfo=TRUE}
  colnames(data)
  key <- data[,c("Sequence", "Modifications")]
  Abu <- data[,grep("Abundance", colnames(data))]
  if (PeakInfo) {
    Peak <- data[,c('RT','m/z')]
  }
  Description <- data[,c('Accession','Description')]
  if (PeakInfo) {
    Peptide<-list(key=key,Abu=Abu,Peak=Peak,
                  Description=Description)
  }
  else {
    Peptide<-list(key=key,Abu=Abu,Description=Description)
  }
  Peptide$Abu[is.na(Peptide$Abu)] <- 0
  return(Peptide)
}
