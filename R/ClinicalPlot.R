#' Plot calibre curve of a lasso cox survival model
#'
#' @param cal rms::calibrate object
#' @param xlim_min min value of x axis
#' @param ylim_max max value of x axis
#' @param ylim_min min value of y axis
#' @param ylim_max max value of y axis
#' @param xlab label of x axis
#' @param ylab label of y axis
#' @return calibre line plot
#' @author Yangming si
#' @examples
#' library(rms)
#' set.seed(1)
#' n <- 200
#' d.time <- rexp(n)
#' x1 <- runif(n)
#' x2 <- factor(sample(c("a", "b", "c"), n, TRUE))
#' f <- cph(Surv(d.time) ~ pol(x1, 2) * x2, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1.5)
#' cal <- calibrate(f, u = 1.5, cmethod = "KM", m = 50, B = 20)
#' calibratePlot(cal)
#' @export
calibratePlot <- function(cal, xlim_min = 0, xlim_max = 1, ylim_min = 0, ylim_max = 1, xlab = "", ylab = "") {
  plot(cal,
    lwd = 2, lty = 1, errbar.col = "black",
    xlim = c(xlim_min, xlim_max), ylim = c(ylim_min, ylim_max),
    xlab = xlab,
    ylab = ylab,
    col = "blue",
    sub = F
  )
  lines(
    x = cal[, "mean.predicted"], y = cal[, "KM"],
    # c('mean.predicted','KM'),
    type = "l", lwd = 3, col = "black",
    pch = 16
  )
  mtext("")
  box(lwd = 0.5)
  abline(0, 1, lty = 3, lwd = 2, col = "black")
}

#' Risk Score Rank Plot of a lasso cox survival model
#'
#' @param riskScore a vector of risk Score
#' @return Risk Score Rank plot
#' @author Yangming si
#' @examples
#' data(riskScore)
#' @export
RiskScoreRankPlot <- function(riskScore) {
  plot(riskScore,
    type = "p",
    pch = 20,
    xlab = "Patients (increasing risk socre)",
    ylab = "Risk score",
    col = c(
      rep("green", lowLength),
      rep("red", highLength)
    )
  )
  abline(
    h = median(riskScore),
    v = lowLength, lty = 2
  )
}

#' Survival curve plot
#'
#' @param fit model from survival::survfit
#' @param geneName gene Name
#' @param pval p Value
#' @return Survival curve plot
#' @author Yangming si
#' @export
SurvivalCurvePlot <- function(fit, geneName, pval = pval) {
  plot(fit,
    col = c("blue", "red"), xlab = "Time (months)", ylab = "Overall survival",
    main = paste(geneName, "(p=", pval, ")", sep = ""), mark.time = T, ylim = c(0, 1.1),
    lwd = 2, cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, font = 1.2
  )
  legend("topright", c(
    paste("Low expression"),
    paste("High expression")
  ),
  col = c("blue", "red"), bty = "n", lwd = 2, cex = 0.8
  )
}
