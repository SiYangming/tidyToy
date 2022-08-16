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
mean2 <- function(x, group = group) {
  tapply(x, group, mean)
}

#' Calculate mean and standard deviation of data frame
#'
#' @param data A data frame
#' @param varname column name of data to calculate
#' @param groupnames column name of data represent different group
#' @return mean and sd of groupname
#' @author Yangming si
#' @examples
#' data("clinical")
#' data_summary(clinical, "futime", "gender")
#' @export

data_summary <- function(data, varname, groupnames) {
  data_sum <- plyr::ddply(data, groupnames,
    .fun = .summary_func,
    all_of(varname)
  )
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

#' Calculate mean and standard deviation
#'
#' @param x A list contain values
#' @param col Name of sub-list
#' @return A vector of mean and standard deviation
#' @author Yangming si
#' @examples
#' See data_summary function
#' @export
.summary_func <- function(x, col) {
  c(
    mean = mean(x[[col]], na.rm = TRUE),
    sd = sd(x[[col]], na.rm = TRUE)
  )
}

#' PValue format
#'
#' @param pValue a vector of p value
#' @return formatted p value
#' @author Yangming si
#' @examples
#' library(qvalue)
#' data(hedenfalk)
#' p <- hedenfalk$p[1]
#' pValueFormat(p)
#' @export

pValueFormat <- function(pValue) {
  pval <- 0
  if (pValue > 0.05) {
    pval <- round(as.numeric(pValue), 3)
  }
  if (pValue < 0.05) {
    pval <- signif(as.numeric(pValue), 4)
    pval <- format(pval, scientific = TRUE)
  }
  return(pval)
}

#' Convert a list to data frame by row or column
#'
#' @param pValue a vector of p value
#' @param by row pr col
#' @return formatted p value
#' @author Yangming si
#' @examples
#' l <- list(rnorm(10), runif(10), rnorm(10, 1))
#' df <- List2DataFrame(l)
#' @export
List2DataFrame <- function(x, by = "row") {
  if (by %in% c("row", "col")) {
    if (by == "row") {
      return(do.call(rbind, x))
    }
    if (by == "col") {
      return(do.call(cbind, x))
    }
  } else {
    print("row or col.")
  }
}
