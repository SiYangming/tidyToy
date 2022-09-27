#' Batch compute Logistic regression model
#'
#' @param pkgs Package names to install
#' @param cran_mirror CRAN Mirror website
#' @param bioc_mirror Bioconductor Mirror website
#' @return NULL
#' @author Yangming si
#' @examples
#' pkgs <- c("tidyverse", "limma", "affy", "oligo", "lumi")
#' pkgs_in(pkgs)
#' @export
pkgs_in <- function(pkgs, cran_mirror = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/", bioc_mirror = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/") {
  # set Chinese mirror
  options("repos" = c(CRAN = cran_mirror))
  options("BioC_mirror" = bioc_mirror)

  # install BiocManager first，skip when installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", ask = F, update = F)
  }

  # 过install stringr first，skip when installed
  if (!requireNamespace("stringr", quietly = TRUE)) {
    install.packages("stringr", ask = F, update = F)
  }

  # remove duplicated，recognize github packages
  pkgs <- unique(pkgs)
  pkgs2 <- pkgs
  logi <- stringr::str_detect(pkgs2, "/")
  pkgs2[logi] <- stringr::str_match(pkgs2[logi], ".*/(.*)$")[, 2]

  # Packages in pkgs not yet installed
  new <- !(sapply(pkgs2, requireNamespace, quietly = T))

  # print pkgs to install
  if (sum(new) > 0) {
    cat("pkgs to install: ", pkgs[new], "\n")
  } else {
    cat("All pkgs already installed \n")
  }

  # install pkgs
  if (any(new)) BiocManager::install(pkgs[new], ask = F, update = F)
}
