#' Batch compute Correlation between two Dataframe
#'
#' @param df1 A data frame
#' @param df2 A data frame
#' @param CorrelationCutoff Correlation Cutoff value
#' @param PvalueCutoff Pvalue Cutoff
#' @param alternative one of "two.sided", "greater" or "less"
#' @param method One of "pearson", "kendall", or "spearman", can be abbreviated
#' @return A tibble of Correlation result
#' @author Yangming si
#' @examples
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' CorTestBatch(df1, df2)
#' @export

CorTestBatch <- function(df1, df2, CorrelationCutoff = 0.5, PvalueCutoff = 0.05, alternative = "two.sided", method = "pearson"){
  # correlation testing
  cor_test <- list()
  for (i in colnames(df1)) {
    for (j in colnames(df2)) {
      x <- as.numeric(df1[, i])
      y <- as.numeric(df2[, j])
      reslut <- list(cor.test(x, y, method = method))
      cor_test <- append(cor_test, reslut)
    }
  }
  # extract correlation test result
  cor_test_lst <- lapply(1:length(cor_test), function(i) {
    broom::tidy(cor_test[[i]])
  })

  cor_test_df <- do.call(bind_rows, cor_test_lst) %>%
    mutate(
      gene1 = rep(colnames(df1), each = ncol(df2)),
      gene2 = rep(colnames(df2), times = ncol(df1)),
      .before = estimate,
      Significance = dplyr::if_else(abs(estimate) >= CorrelationCutoff &
                                      p.value < PvalueCutoff,
                                    dplyr::if_else(estimate > 0, "Postive", "Negative"),
                                    "Not Sig")
    )
  return(cor_test_df)
}
