#' Box plot for immunity Results
#'
#' @param genename gene Name
#' @param cell Cell Name: One or More
#' @param gene_mat Gene Expression Matrix
#' @param cell_matrix Cell signal Matrix
#' @return A ggplot2 object of box plot
#' @author Yangming si
#' @examples
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' ggboxplot2("gene1", "gene3", df1, df2)
#' @export

ggboxplot2 <- function(genename, cell, gene_mat, cell_matrix) {
  gene <- as.numeric(gene_mat[, genename])
  group <- if_else(gene > median(gene),
                   paste0(genename, "-High"),
                   paste0(genename, "-Low"))
  group <- factor(group, levels = c(paste0(genename, "-Low"),
                                    paste0(genename, "-High")))
  if (length(cell) > 1) {
    box <- reshape2::melt(data.frame(group = group,
                           cell_matrix[,cell]))
  } else {
    box <- reshape2::melt(data.frame(group = group,
                           cell_matrix[,cell]))
    box$variable <- cell
  }

  box$variable <- str_replace_all(box$variable, "\\.", " ")
  p <- ggpubr::ggboxplot(box,
                 x = "variable", y = "value", color = "group",
                 palette = c("#00AFBB", "#E7B800")
  ) +
    # easy_rotate_x_labels(angle = 45, teach = TRUE)
    # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggpubr::stat_compare_means(aes(group = group), label = "p.signif") +
    xlab("") + ylab("")
  return(p)
}

#' Scatter plot for immunity Results
#'
#' @param genename gene Name
#' @param cell Cell Name
#' @param gene_mat Gene Expression Matrix
#' @param cell_matrix Cell signal Matrix
#' @param trend pos or neg
#' @return A ggplot2 object of scatter plot
#' @author Yangming si
#' @examples
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' ggscatter2("gene1", "gene3", df1, df2)
#' @export

ggscatter2 <- function(genename, cell, gene_mat, cell_matrix, trend = "pos") {
  x <- gene_mat[, genename]
  scatter <- data.frame(x = x, y = cell_matrix[,cell])
  if (trend == "pos") {
    label.x = min(scatter$x) + 0.5
    label.y = max(scatter$y) - 0.03
  }
  if (trend == "neg") {
    label.x = max(scatter$x) - 0.5
    label.y = max(scatter$y) - 0.03
  }
  p <- ggpubr::ggscatter(scatter,
                 x = "x", y = "y",
                 fill = "black",
                 color = "black", shape = 21, size = 1, # Points color, shape and size
                 # add = "reg.line", # Add regressin line
                 xlab = genename, ylab = cell,
                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg.line
                 conf.int = TRUE, # Add confidence interval
                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                 cor.coeff.args = list(method = "pearson",
                                       label.x = label.x,
                                       label.y = label.y,
                                       label.sep = "\n")
  ) +
    ggplot2::geom_smooth(method = "lm")
  return(p)
}
