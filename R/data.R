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

#' immunity cell signature for ssGSEA
#'
#' A dataset containing gene signature in 24 immunity cell
#'
#' @format A list
#' \describe{
#'   \item{aDC}{gene signature in aDC cell}
#'   \item{B cells}{gene signature in B cells cell}
#'   \item{CD8 T cells}{gene signature in CD8 T cells cell}
#'   \item{Cytotoxic cells}{gene signature in Cytotoxic cells cell}
#'   \item{DC}{gene signature in DC cell}
#'   \item{Eosinophils}{gene signature in Eosinophils cell}
#'   \item{iDC}{gene signature in iDC cell}
#'   \item{Macrophages}{gene signature in Macrophages cell}
#'   \item{Mast cells}{gene signature in Mast cells cell}
#'   \item{Neutrophils}{gene signature in Neutrophils cell}
#'   \item{NK CD56bright cells}{gene signature in NK CD56bright cells cell}
#'   \item{NK CD56dim cells}{gene signature in NK CD56dim cells cell}
#'   \item{NK cells}{gene signature in NK cells cell}
#'   \item{pDC}{gene signature in pDC cell}
#'   \item{T cells}{gene signature in T cells cell}
#'   \item{T helper cells}{gene signature in T helper cells cell}
#'   \item{Tcm}{gene signature in Tcm cell}
#'   \item{Tem}{gene signature in Tem cell}
#'   \item{TFH}{gene signature in TFH cell}
#'   \item{Tgd}{gene signature in Tgd cell}
#'   \item{Th1 cells}{gene signature in Th1 cells cell}
#'   \item{Th17 cells}{gene signature in Th17 cells cell}
#'   \item{Th2 cells}{gene signature in Th2 cells cell}
#'   \item{TReg}{gene signature in TReg cell}
#' }
#' @name immunity24
#' @source table S1 - https://doi.org/10.1016/j.immuni.2013.10.003
#' @docType data
#' @keywords datasets
#' @usage data(immunity24)
NULL
