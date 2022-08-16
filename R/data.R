#' riskScore data for plot
#'
#' A dataset containing Cox model predicted results
#'
#' @format A data frame
#' \describe{
#'   \item{id}{TCGA patient Identifier}
#'   \item{futime}{survival time}
#'   \item{fustat}{survival status}
#'   \item{age}{How old is patient}
#'   \item{grade}{cancer grade}
#'   \item{stage}{tumor stage}
#'   \item{T}{tumor T stage}
#'   \item{M}{tumor M stage}
#'   \item{N}{tumor N stage}
#'   \item{riskScore}{Cox model predicted result}
#' }
#' @name clinical
#' @source TCGA-STAD
#' @docType data
#' @keywords datasets
#' @usage data(clinical)
NULL
