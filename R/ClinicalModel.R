#' survival cox model
#'
#' @param x a data frame of survival data
#' @return A data frame of survival model
#' @author Yangming si
#' @examples
#' data(clinical)
#' x <- as.data.frame(clinical)[-1][1:3]
#' SurvivalModel(x)
#' @export
SurvivalModel <- function(x) {
  cox <- coxph(Surv(futime, fustat) ~ ., data = x)
  coxSummary <- summary(cox)
  coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
  gene <- as.numeric(x[, 3, drop = TRUE])
  med <- median(gene)
  if (med != 0) {
    rt1 <- filter(x, gene >= med)
    rt2 <- filter(x, gene < med)
    surTab1 <- summary(survfit(Surv(futime, fustat) ~ 1, data = rt1))
    surTab2 <- summary(survfit(Surv(futime, fustat) ~ 1, data = rt2))
    medianTab1 <- surTab1$table
    medianTab2 <- surTab2$table
    diff <- survdiff(Surv(futime, fustat) ~ (gene > med), data = x)
    fit <- survfit(Surv(futime, fustat) ~ (gene > med), data = x)
    pValue <- 1 - pchisq(diff$chisq, df = 1)
    coef <- as.vector(coxSummary$coefficients)
    names(coef) <- colnames(coxSummary$coefficients)
    conf <- as.vector(coxSummary$conf.int)
    names(conf) <- colnames(coxSummary$conf.int)
    res <- c(
      gene = colnames(x)[3],
      coef,
      conf,
      KM = pValue,
      H_med = medianTab1["median"],
      H_0.95LCL = medianTab1["0.95LCL"],
      H_0.95UCL = medianTab1["0.95UCL"],
      L_med = medianTab2["median"],
      L_0.95LCL = medianTab2["0.95LCL"],
      L_0.95UCL = medianTab2["0.95UCL"]
    )
  }
  return(res)
}
#' unicox multicox model
#'
#' @param mydata a data frame of survival data
#' @return cox result table
#' @author Yangming si
#' @examples
#' data(clinical)
#' x <- as.data.frame(clinical)[-1][1:3]
#' cox8000(x)
#' @export
cox8000 <- function(mydata0) {
  library(survival)
  result0 <- list()
  HR0 <- list()
  HRCOEF0 <- list()
  for (i in 3:ncol(mydata0))
  {
    result0[[i - 2]] <- anova(coxph(Surv(futime, fustat) ~ mydata0[, i], data = mydata0, ties = "breslow"))$Pr[2]
    HR0[[i - 2]] <- summary(coxph(Surv(futime, fustat) ~ mydata0[, i], data = mydata0, ties = "breslow"))$conf.int[1, ]
    HRCOEF0[[i - 2]] <- summary(coxph(Surv(futime, fustat) ~ mydata0[, i], data = mydata0, ties = "breslow"))$coefficients[1, ]
  }
  data0 <- t(mydata0[-c(1:2)])
  pval0 <- unlist(result0)
  gene0 <- rownames(data0)
  HRR0 <- do.call(rbind, lapply(HR0, `[`, c(1:4))) ####### 非常重要
  HRR0 <- data.frame(HRR0)
  names(HRR0) <- c("exp(coef)", "exp(-coef)", "lower .95", "upper .95") ######
  HRCOEFF0 <- do.call(rbind, lapply(HRCOEF0, `[`, c(1:5))) ####### 非常重要
  HRCOEFF0 <- data.frame(HRCOEFF0)
  names(HRCOEFF0) <- c("coef1", "exp(coef)1", "se(coef)1", "z1", "Pvalue1") # coef coef   coef z Pr(>|z|)
  resultsata0 <- data.frame(gene0, pval0, HRR0, HRCOEFF0)
  return(resultsata0)
}
