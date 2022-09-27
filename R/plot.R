#' Box plot for One vs More Variables: immunity Results (p starred)
#'
#' @param col1 column names of input data frame, only one
#' @param col2 column names of input data frame, One or More
#' @param mat Expression Matrix or data frame, not tibble
#' @param group group column name or null, default NULL
#' @return A ggplot2 object of box plot
#' @author Yangming si
#' @examples
#' # example 1: Low High group
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' df <- cbind(df1, df2)
#' ggboxplot2("gene1", "gene3", df)
#' # example 2: multi group
#' df1 <- data.frame(gene1 = rnorm(30), gene2 = rnorm(30), gene3 = rnorm(30))
#' df2 <- data.frame(group = rep(c(1, 2, 3), each = 10))
#' df <- cbind(df1, df2)
#' ggboxplot2("gene1", "group", df, group = "group")
#' @export

ggboxplot2 <- function(col1, col2, mat, group = NULL) {
  gene <- as.numeric(mat[, col1])
  if (is.null(group)) {
    group <- if_else(gene > median(gene),
      paste0(col1, "-High"),
      paste0(col1, "-Low")
    )
    group <- factor(group, levels = c(
      paste0(col1, "-Low"),
      paste0(col1, "-High")
    ))
    if (length(col2) > 1) {
      box <- reshape2::melt(data.frame(
        group = group,
        mat[, col2]
      ))
    } else {
      box <- reshape2::melt(data.frame(
        group = group,
        mat[, col2]
      ))
      box$variable <- col2
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
  } else {
    box <- mat[, c(col1, group)]
    my_comparisons <- Group2Comparison(unique(sort(as.character(box[, group]))))
    p <- ggpubr::ggboxplot(box,
      x = group,
      y = col1, color = group
    ) +
      ggpubr::stat_compare_means(ggplot2::aes(group = group),
        comparisons = my_comparisons,
        label = "p.signif"
      ) +
      ggplot2::xlab(group) + ggplot2::ylab(col1)
  }


  return(p)
}

#' Box plot for One vs More Variables: immunity Results (pvalue)
#'
#' @param col1 column names of input data frame, only one
#' @param col2 column names of input data frame, One or More
#' @param mat Expression Matrix or data frame, not tibble
#' @param group group column name or null, default NULL
#' @return A ggplot2 object of box plot
#' @author Yangming si
#' @examples
#' # example 1: Low High group
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' df <- cbind(df1, df2)
#' ggboxplot2("gene1", "gene3", df)
#' # example 2: multi group
#' df1 <- data.frame(gene1 = rnorm(30), gene2 = rnorm(30), gene3 = rnorm(30))
#' df2 <- data.frame(group = rep(c(1, 2, 3), each = 10))
#' df <- cbind(df1, df2)
#' ggboxplot2("gene1", "group", df, group = "group")
#' @export

ggboxplot3 <- function(col1, col2, mat, group = NULL) {
  gene <- as.numeric(mat[, col1])
  if (is.null(group)) {
    group <- if_else(gene > median(gene),
                     paste0(col1, "-High"),
                     paste0(col1, "-Low")
    )
    group <- factor(group, levels = c(
      paste0(col1, "-Low"),
      paste0(col1, "-High")
    ))
    if (length(col2) > 1) {
      box <- reshape2::melt(data.frame(
        group = group,
        mat[, col2]
      ))
    } else {
      box <- reshape2::melt(data.frame(
        group = group,
        mat[, col2]
      ))
      box$variable <- col2
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
  } else {
    box <- mat[, c(col1, group)]
    my_comparisons <- Group2Comparison(unique(sort(as.character(box[, group]))))
    p <- ggpubr::ggboxplot(box,
                           x = group,
                           y = col1, color = group
    ) +
      ggpubr::stat_compare_means(ggplot2::aes(group = group),
                                 comparisons = my_comparisons,
                                 label = "p.format"
      ) +
      ggplot2::xlab(group) + ggplot2::ylab(col1)
  }


  return(p)
}

#' Covert Group vector to a comparison list
#'
#' @param group group character vector
#' @return A group comparison list
#' @author Yangming si
#' @examples
#' group <- as.character(1:3)
#' Group2Comparison(group)
#' @export

Group2Comparison <- function(group) {
  return(combn(group, 2, simplify = FALSE))
}

#' Scatter plot for Two Variables
#'
#' @param col1 column names of input data frame
#' @param col2 Ccolumn names of input data frame
#' @param mat Expression Matrix or data frame, not tibble
#' @param trend pos or neg
#' @param method one of "pearson", "kendall", or "spearman"
#' @return A ggplot2 object of scatter plot
#' @author Yangming si
#' @examples
#' df1 <- data.frame(gene1 = rnorm(50), gene2 = rnorm(50))
#' df2 <- data.frame(gene3 = rnorm(50), gene4 = rnorm(50))
#' df <- cbind(df1, df2)
#' ggscatter2("gene1", "gene3", df)
#' @export

ggscatter2 <- function(col1, col2, mat, trend = "pos", method = "pearson") {
  scatter <- data.frame(x = mat[, col1], y = mat[, col2])

  trend <- ifelse(cor.test(scatter$x, scatter$y,
    method = method
  )$estimate > 0,
  "pos", "neg"
  )
  if (trend == "pos") {
    label.x <- mean(scatter$x)
    label.y <- max(scatter$y) - mean(scatter$y) / 2
  }
  if (trend == "neg") {
    label.x <- mean(scatter$x)
    label.y <- max(scatter$y) - mean(scatter$y) / 2
  }
  p <- ggpubr::ggscatter(scatter,
    x = "x", y = "y",
    fill = "black",
    color = "black", shape = 21, size = 1, # Points color, shape and size
    # add = "reg.line", # Add regressin line
    xlab = col1, ylab = col2,
    add.params = list(color = "blue", fill = "lightgray"), # Customize reg.line
    conf.int = TRUE, # Add confidence interval
    cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
    cor.coeff.args = list(
      method = method,
      label.x = label.x,
      label.y = label.y,
      label.sep = "\n"
    )
  ) +
    ggplot2::geom_smooth(method = "lm")
  return(p)
}

#' Plots wordcloud.
#'
#' @description Plots word cloud of mutated genes or altered cytobands with size proportional to the event frequency.https://gist.github.com/PoisonAlien/3f8752a89c1d63f64afe55b441621223
#' @param input an \code{\link{MAF}} or \code{\link{GISTIC}} object generated by \code{\link{read.maf}} or \code{\link{readGistic}}
#' @param minMut Minimum number of samples in which a gene is required to be mutated.
#' @param col vector of colors to choose from.
#' @param top Just plot these top n number of mutated genes.
#' @param genesToIgnore Ignore these genes.
#' @param ... Other options passed to \code{\link{wordcloud}}
#' @return nothing.
#' @examples
#' library(maftools)
#' laml.input <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read.maf(maf = laml.input, useAll = FALSE)
#' geneCloud(input = laml, minMut = 5, top = 5)
#' @importFrom wordcloud wordcloud
#' @export


geneCloud <- function(input, minMut = 3, col = NULL, top = NULL, genesToIgnore = NULL, ...) {
  if (class(input)[1] == "GISTIC") {
    gs <- getCytobandSummary(input)

    col <- ifelse(test = gs$Variant_Classification %in% "Amp", yes = "red", no = "blue")
    wordcloud::wordcloud(
      words = gs[, Cytoband], freq = -log10(gs[, qvalues]),
      colors = col, ordered.colors = TRUE, rot.per = 0.2, ...
    )
  } else {
    gs <- getGeneSummary(x = input)

    # Either choose top n genes or > minMut
    if (!is.null(top)) {
      gs <- gs[1:top]
    } else {
      gs <- gs[gs$MutatedSamples >= minMut, ]
    }

    # If any genes to ignore
    if (!is.null(genesToIgnore)) {
      gs <- gs[!Hugo_Symbol %in% genesToIgnore]
    }

    if (is.null(col)) {
      col <- c(RColorBrewer::brewer.pal(12, name = "Paired"), RColorBrewer::brewer.pal(9, name = "Set1"), "black")
    }

    wordcloud::wordcloud(
      words = gs[["Hugo_Symbol"]], freq = gs[["MutatedSamples"]],
      colors = col, random.color = TRUE, rot.per = 0.2, ...
    )
  }
}

#' Plot Volcano.
#'
#' @description Volcano Plot for DEGs
#' @param DEGs generated by \code{\link{GetArrayDEGs}}
#' @param logfc Log Fold Change
#' @param pvalue P value
#' @param col vector of colors to choose from.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param v Horizontalism line
#' @param h Vertical line
#' @param line_col line color
#' @param ... Other options passed to \code{\link{plot}}
#' @return nothing.
#' @examples
#' library(CLL)
#' data("CLLbatch")
#' mat <- exprs(CLLbatch)
#' mat <- LOG2Transform(mat)
#' group <- rep(c("Control", "Treat"), each = ncol(mat) / 2)
#' design <- GetDesign(group, mat)
#' contrasts <- "Treat-Control"
#' deg <- GetArrayDEGs(mat, contrasts, design, logFCcutoff = 1, PvalueCutoff = 0.05)
#' col <- ifelse(deg$P.Value <= 0.05 & abs(deg$logFC) > 1, "blue", "gray")
#' VolcanoPlot(deg, col = col, main = "Volcano Plot")
#' colour1 <- ifelse(allDiff$P.Value <= 0.05 & abs(allDiff$logFC) > 1, "green", "gray")
#' colour2 <- ifelse(allDiff$P.Value > 0.05, "gray", ifelse(allDiff$logFC > 1, "red", ifelse(allDiff$logFC < -1, "green", "gray")))
#' @export

VolcanoPlot <- function(DEGs, logfc = "logFC", pvalue = "P.Value",
                        col = NULL, xlab = "log2 fold chage",
                        ylab = "-log10(p-value)",
                        v = c(-1, 1), h = -log10(0.05),
                        line_col = "green", ...) {
  # par(mar(3,3,2,1), mgp = c(1.6,0.6,0), tck = .01)
  xMax <- max(abs(DEGs[[logfc]]))
  yMax <- max(-log10(-log10(DEGs[[pvalue]])))
  if (is.null(col)) {
    col <- "grey"
  }
  plot(DEGs[[logfc]], -log10(DEGs[[pvalue]]),
    pch = 16, # 设置点的类型
    cex = 0.5, # 设置点的大小
    cex.lab = 1.2, # 设置坐标轴名称的字体大小
    font.lab = 2, # 设置坐标轴名称的字体类型，粗体
    col = col, # 设置点的颜色
    xlim = c(-xMax, xMax), ylim = c(0, yMax), # Set limits
    xlab = xlab, # 设置坐标轴标签
    ylab = ylab, ...
  )
  abline(
    v = c(-1, 1),
    h = -log10(0.05), col = line_col,
    lty = 2, lwd = 3
  )
  diffSub <- subset(DEGs, P.Value < 0.05 & logFC > 1)
  points(diffSub$logFC,
         -log10(diffSub$P.Value),
    pch = 20, col = "red", cex = 0.4
  )
  diffSub <- subset(DEGs, P.Value < 0.05 & logFC < -1)
  points(diffSub$logFC,
         -log10(diffSub$P.Value),
    pch = 20, col = "blue", cex = 0.4
  )
  # 在图上标记基因
  ids <- order(deg[["logFC"]], decreasing = TRUE)
  # grep("CD36", DiffGene$symbol)
  # x[grep("CD36", DiffGene$symbol)]
  # y[grep("CD36", DiffGene$symbol)]
  text(DEGs[[logfc]][ids[1]], -log10(DEGs[[pvalue]][ids[1]]), rownames(DEGs)[ids[1]])
  text(DEGs[[logfc]][ids[length(ids)]], -log10(DEGs[[pvalue]][ids[length(ids)]]), rownames(DEGs)[ids[length(ids)]])
}

#' Plot PCA.
#'
#' @description PCA Plot for Gene Expression
#' @param expr Gene Expression Matrix
#' @param ntop How many genes to use
#' @param group group vector
#' @param show_name show name or not.
#' @return ggplot object.
#' @examples
#' library(CLL)
#' data("CLLbatch")
#' mat <- exprs(CLLbatch)
#' group <- rep(c("Control", "Treat"), each = ncol(mat) / 2)
#' PCA_new(mat, group = group)
#' @export
PCA_new <- function(expr, ntop = 500, group, show_name = F){
  library(ggplot2)
  library(ggrepel)
  object <- expr
  rv <- genefilter::rowVars(object)
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(object[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[, 1],
                  PC2 = pca$x[, 2],
                  group = group,
                  name = colnames(object))
  attr(d, "percentVar") <- percentVar[1:2]
  if (show_name) {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
      geom_text_repel(aes(label = name),
                      size = 3,
                      segment.color = "black",
                      show.legend = FALSE )
  } else {
    ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = "group")) +
      geom_point(size = 2) +
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
  }
}
