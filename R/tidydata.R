#' filter CIBERSORT Result
#'
#' @param results CIBERSORT Result matrix
#' @param percent remain A column in a vector containing a proportion of zero less than 'percent', default 0.5
#' @param CorrelationCutoff Correlation Cutoff value
#' @param PvalueCutoff Pvalue Cutoff
#' @param alternative one of "two.sided", "greater" or "less"
#' @param method One of "pearson", "kendall", or "spearman", can be abbreviated
#' @return filtered matrix
#' @author Yangming si
#' @examples
#' library(CIBERSORT)
#' data(LM22)
#' data(mixed_expr)
#' results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
#' filtered_results <- filter_CIBERSORT(results)
#' @export

filter_CIBERSORT <- function(results, percent = 0.5){
  results <- results[, -c((ncol(results) - 2):ncol(results))]

  remain_cells <- sapply(colnames(results), function(x){
    sum(near(results[, x], 0)) / nrow(results) < percent
  })
  return(results[, remain_cells])
}
