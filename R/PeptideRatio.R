#' Title PeptideRatio
#' This function calculates pairwise ratio of peptides between two subpopulations
#'
#' @param Peptide Quantitative peptide information without untrusted data and contaminant data, and without outlier samples
#' @param numerator the type of cells that play the numerator when calculating the peptide ratio
#'
#' @return A data list of peptide information with a table of peptide DE-scores
#' @examples
#' Peptide<-Initiation(test_raw)
#' Peptide<-FilterData(Peptide,type=test_type,'A','B')
#' Peptide<-FilterSample(Peptide)
#' Peptide<-Normalization(Peptide)
#' Peptide<-PeptideRatio(Peptide)
PeptideRatio<-function (Peptide, numerator = Peptide$type[1])
{
  A <- numerator
  B <- unique(Peptide$type)[unique(Peptide$type) != A]
  PairRatio <- function(x, y) {
    if (max(y) == 0) {
      return(10)
    }
    else {
      R <- c()
      for (i in seq_len(length(x))) {
        if (x[i] == 0) {
          R <- c(R, sapply(y, function(a) if (a == 0) (-1) else (x[i]/a)))
        }
        else {
          R <- c(R, sapply(y, function(a) if (a == 0) (Inf) else (x[i]/a)))
        }
      }
      R[R == Inf] <- max(stats::na.omit(as.numeric(unlist(R[R !=
                                                              Inf]))))
      return(if (max(x) > 0) {
        stats::median(stats::na.omit(R[R>=0]))
      } else {
        0
      })
    }
  }
  Ratio <- apply(Peptide$Abu, 1, function(x) log2(PairRatio(x[Peptide$type ==
                                                                A], x[Peptide$type == B])))
  Ratio[Ratio < (-2)] <- (-2)
  Ratio[Ratio > 2] <- 2
  CalculateP <- function(x, y) {
    if (max(length(x), length(y)) == 3) {
      return(stats::t.test(x, y, var.equal = F, warning = F)$p.value)
    }
    else {
      options(warn = -1)
      return(stats::wilcox.test(x, y)$p.value)
    }
  }
  Pval <- apply(Peptide$Abu, 1, function(x) CalculateP(x[Peptide$type ==
                                                           A], x[Peptide$type == B]))
  PeptideRatio <- data.frame(row.names = rownames(Peptide$Abu),
                             Ratio = Ratio, Pval = Pval)
  Peptide <- list(key = Peptide$key, Abu = Peptide$Abu[, Peptide$type ==
                                                         A | Peptide$type == B], Description = Peptide$Description,
                  type = Peptide$type, PeptideRatio = PeptideRatio)
  return(Peptide)
}
