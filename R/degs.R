#' Convert expression data frame to no duplicated matrix
#'
#' @param df A data frame
#' @param nodups logistic value, TRUE value will remove duplicated row in df
#' @return A matrix
#' @author Yangming si
#' @examples
#' df <- data.frame(gene = c("gene1", "gene2", "gene1", "gene2"),
#'                  sample1 = c(1,2,3, 5), sample2 = c(11,12,13, 6),
#'                  sample1 = c(7, 8, 9, 9), sample2 = c(4, 5, 6, 10))
#' DF2matrix(df)
#' @export
DF2matrix <- function(df, nodups = TRUE)
{
  rt <- as.matrix(df)
  rownames(rt) <- rt[,1]
  exp <- rt[,2:ncol(rt)]
  dimnames <- list(rownames(exp), colnames(exp))
  data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
  if (nodups) {
    return(limma::avereps(data))
  }else{
    return(data)
  }
}

#' Convert group vector to design matrix
#'
#' @param group character vector of group
#' @return A group design matrix
#' @author Yangming si
#' @examples
#' df <- data.frame(gene = c("gene1", "gene2", "gene1", "gene2"),
#'                  sample1 = c(1,2,3, 5), sample2 = c(11,12,13, 6),
#'                  sample1 = c(7, 8, 9, 9), sample2 = c(4, 5, 6, 10))
#' mat <- DF2matrix(df)
#' group <- rep(c("Control", "Treat"), each = 2)
#' GetDesign(group, mat)
#' @export
GetDesign <- function(group, mat)
{
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(mat)
  return(design)
}

#' Automatic determination of whether to log transform the array expression matrix
#'
#' @param mat Gene expression matrix
#' @return A raw or log-transformed(if needed) Data Frame
#' @author Yangming si
#' @examples
#' library(CLL)
#' data("CLLbatch")
#' mat <- exprs(CLLbatch)
#' mat1 <- LOG2Transform(mat)
#' @export
LOG2Transform <- function(mat){
  # log2 transform
  qx <- as.numeric(quantile(mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (LogC) {
    mat[which(mat <= 0)] <- NaN
    mat <- log2(mat)
  }
  return(mat)
}


#' Differential expression analysis on microarrays using limma
#'
#' @param mat Gene expression matrix
#' @param contrasts Set comparison group, e.g. treat group - control group
#' @param design group design matrix from `model.matrix()` function
#' @param logFCcutoff cutoff value of logFC
#' @param PvalueCutoff cutoff value of p
#' @param group "Two" or "Multi"
#' @return A Data Frame of DEGs result
#' @author Yangming si
#' @examples
#' library(CLL)
#' data("CLLbatch")
#' mat <- exprs(CLLbatch)
#' mat <- LOG2Transform(mat)
#' group <- rep(c("Control", "Treat"), each = ncol(mat)/2)
#' design <- GetDesign(group, mat)
#' contrasts <- "Treat-Control"
#' deg  <- GetArrayDEGs(mat, contrasts, design, logFCcutoff = 1, PvalueCutoff = 0.05)
#' head(deg)
#' @export
GetArrayDEGs <- function(mat, contrasts, design, logFCcutoff = 1, PvalueCutoff = 0.05, group = "Two")
{
  # Average of rows with the same gene name
  mat <- limma::avereps(mat)
  # Set Group for Differential expression analysis
  cont.matrix <- limma::makeContrasts(contrasts = contrasts,
                                      levels = design)
  # Calculating model covariance
  fit <- limma::lmFit(mat, design) # fit linear model
  fit <- limma::contrasts.fit(fit, cont.matrix)
  fit <- limma::eBayes(fit)
  DEGs <- limma::topTable(fit, adjust = "fdr", n = Inf)
  DEGs[["Significance"]] <- ifelse(DEGs$logFC >= logFCcutoff & DEGs$P.Value < PvalueCutoff, "Up",
                                   ifelse(DEGs$logFC <= -logFCcutoff & DEGs$P.Value < PvalueCutoff, "Down",
                                          "Not Sig"))
  # DEGs <- DEGs %>%
  #   mutate(Significance = if_else(logFC >= 1 & P.Value < 0.05, "Up",
  #                                 if_else(logFC <= -1 & P.Value < 0.05, "Down",
  #                                         "Not Sig")
  #   ))
  return(DEGs)
}
