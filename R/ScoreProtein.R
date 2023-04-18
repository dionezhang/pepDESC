#' Title ScoreProtein
#' This function calculates the protein DE-scores based on peptide DE-scores and peptide abundances
#'
#' @param pvalcutoff Maximum p-value for hypothesis test for differential expressed peptides.
#' @param Peptide Quantitative peptide information with a table of peptide DE-scores
#'
#' @return a table of protein DE-score
#' @examples
#' Peptide<-Initiation(test_raw)
#' Peptide<-FilterData(Peptide,type=test_type,'A','B')
#' Peptide<-FilterSample(Peptide)
#' Peptide<-Normalization(Peptide)
#' Peptide<-PeptideRatio(Peptide)
#' Result<-ScoreProtein(Peptide)
ScoreProtein<-function (Peptide, pvalcutoff = 0.05)
{
  Protein <- Peptide$Description$Accession
  Protein <- unlist(strsplit(Protein, "; "))
  Protein <- Protein[!duplicated(Protein)]
  CorScore <- function(pep, Abu) {
    if (dim(pep)[1] == 1) {
      return(2)
    }
    else {
      abu <- Abu[rownames(pep), ]
      temp <- stats::cor(t(abu))
      temp[temp < 0] <- 0
      abu <- apply(abu, 1, mean)
      abu <- abu * temp
      abu <- apply(abu, 2, sum)
      return(abu)
    }
  }
  result <- c()
  pb<-progress::progress_bar$new(total = length(Protein), clear = FALSE)
  for (i in seq_len(length(Protein))) {
    pb$tick()
    Accession <- Protein[i]
    pep <- rownames(Peptide$Description)[grep(Accession,
                                              Peptide$Description$Accession)]
    pep <- Peptide$PeptideRatio[pep, ]
    pep[pep$Ratio > 1.5, "Ratio"] <- 1.5
    pep[pep$Ratio < (-1.5), "Ratio"] <- (-1.5)
    pep$weight <- CorScore(pep, Peptide$Abu)
    pep$weight <- pep$weight/sum(pep$weight)
    pep <- pep[pep$Pval <= pvalcutoff, ]
    score <- pep$Ratio * pep$weight
    result <- c(result, sum(score))
  }
  result <- data.frame(row.names = Protein, score = result)
  return(result)
}

