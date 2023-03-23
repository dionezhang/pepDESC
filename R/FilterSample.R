#' @title FilterSample
#' @document This function removes outlier samples according to the identification and total abundance.
#' @param Peptide Quantitative peptide processed without untrusted data or contaminant data
#'
#' @return A data list of peptide information without outlier samples
#' @examples
#' Peptide<-Initiation(test_raw)
#' Peptide<-FilterData(Peptide,type=test_type,'A','B')
#' Peptide<-FilterSample(Peptide)
FilterSample<-function (Peptide)
{
  data <- Peptide$Abu
  sum <- apply(data, 2, sum)
  count <- apply(data, 2, function(x) sum(x > 0))
  if (stats::shapiro.test(sum)$p.value > 0.05) {
    sum <- (1/sum)
  }
  if (stats::shapiro.test(count)$p.value > 0.05) {
    count <- (1/count)
  }
  fit_sum <- MASS::fitdistr(sum, "normal")
  fit_count <- MASS::fitdistr(count, "normal")
  rmSample <- names(count[count < (fit_count$estimate[1] -
                                     2 * fit_count$estimate[2])])
  rmSample <- c(rmSample, names(sum[sum < (fit_sum$estimate[1] -
                                             2 * fit_sum$estimate[2])]))
  rmSample <- unique(rmSample)
  print(paste("removed sample =", as.character(rmSample),
              sep = "\r"))
  Peptide$type <- Peptide$type[colnames(Peptide$Abu) %in%
                                 rmSample == F]
  Peptide$Abu <- Peptide$Abu[, colnames(Peptide$Abu) %in%
                               rmSample == F]
  return(Peptide)
}

