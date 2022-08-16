#' Calculate average of different group
#'
#' @param x A data frame or tibble
#' @param group A factor character vector
#' @return A vector of different group average
#' @author Yangming si
#' @examples
#' x <- c(rnorm(10), runif(10), rnorm(10, 1))
#' f <- gl(3, 10)
#' f
#' mean2(x, f)
#' @export
mean2 <- function(x, group = group){
  tapply(x, group, mean)
}
