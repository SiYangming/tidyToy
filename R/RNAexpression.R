#' Convert count matrix to tpm expression
#'
#' @param counts A count matrix from NGS
#' @param effLen effective length from NGS
#' @return A TPM expression matrix
#' @author Yangming si
#' @examples
#' cnts <- c(4250, 3300, 200, 1750, 50, 0)
#' lens <- c(900, 1020, 2000, 770, 3000, 1777)
#' countDf <- data.frame(count = cnts, length = lens)
#' # assume a mean(FLD) = 203.7
#' countDf$effLength <- countDf$length - 203.7 + 1
#' countDf$tpm <- with(countDf, countToTpm(count, effLength))
#' countDf$tpm
#' @export

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

#' Convert count matrix to Fpkm expression
#'
#' @param counts A count matrix from NGS
#' @param effLen effective length from NGS
#' @return A Fpkm expression matrix
#' @author Yangming si
#' @examples
#' cnts <- c(4250, 3300, 200, 1750, 50, 0)
#' lens <- c(900, 1020, 2000, 770, 3000, 1777)
#' countDf <- data.frame(count = cnts, length = lens)
#' # assume a mean(FLD) = 203.7
#' countDf$effLength <- countDf$length - 203.7 + 1
#' countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
#' countDf$fpkm
#' @export

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

#' Convert Fpkm matrix to Tpm expression
#'
#' @param fpkm Fpkm count matrix from NGS
#' @return A Tpm expression matrix
#' @author Yangming si
#' @examples
#' cnts <- c(4250, 3300, 200, 1750, 50, 0)
#' lens <- c(900, 1020, 2000, 770, 3000, 1777)
#' countDf <- data.frame(count = cnts, length = lens)
#' # assume a mean(FLD) = 203.7
#' countDf$effLength <- countDf$length - 203.7 + 1
#' countDf$fpkm <- with(countDf, countToFpkm(count, effLength))
#' with(countDf, all.equal(tpm, fpkmToTpm(fpkm)))
#' @export

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#' Convert count matrix to effective Counts
#'
#' @param counts A count matrix from NGS
#' @param effLen effective length from NGS
#' @return A effective Counts matrix
#' @author Yangming si
#' @examples
#' cnts <- c(4250, 3300, 200, 1750, 50, 0)
#' lens <- c(900, 1020, 2000, 770, 3000, 1777)
#' countDf <- data.frame(count = cnts, length = lens)
#' # assume a mean(FLD) = 203.7
#' countDf$effLength <- countDf$length - 203.7 + 1
#' countDf$effCounts <- with(countDf, countToEffCounts(count, length, effLength))
#' @export

countToEffCounts <- function(counts, len, effLen)
{
  counts * (len / effLen)
}
