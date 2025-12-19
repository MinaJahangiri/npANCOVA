#' npANCOVA: Nonparametric ANCOVA Methods
#'
#' @description
#' The npANCOVA package provides a collection of functions to implement various
#' nonparametric Analysis of Covariance (ANCOVA) methods.
#'
#' @details
#' This package is designed for researchers and statisticians who need to
#' perform ANCOVA when the assumptions of parametric ANCOVA, such as
#' normality of residuals or the homogeneity of regression slopes, are not met.
#'
#' @section Key Functions:
#' The main functions provided by the package are:
#' \describe{
#'   \item{\code{\link{Burnett_Barr}}}{Implements the Burnett and Barr rank-based method.}
#'   \item{\code{\link{Harwell_Serlin}}}{Performs the Harwell and Serlin method using ranked variables.}
#'   \item{\code{\link{Hettmansperger_McKean}}}{Applies rank-based residual analysis.}
#'   \item{\code{\link{McSweeny_Porter}}}{Performs rank-based ANCOVA with and without interaction.}
#'   \item{\code{\link{Puri_Sen_OU}}}{Performs the Puri and Sen method using one covariate with an ubiased variance-covariance matrix.}
#'   \item{\code{\link{Puri_Sen_OB}}}{Performs the Puri and Sen method using one covariate with a biased variance-covariance matrix.}
#'   \item{\code{\link{Puri_Sen_MU}}}{Performs the Puri and Sen method using multiple covariates with an unbiased variance-covariance matrix.}
#'   \item{\code{\link{Puri_Sen_MB}}}{Performs the Puri and Sen method using multiple covariates with a biased variance-covariance matrix.}
#'   \item{\code{\link{Quade}}}{Performs Quade's ANCOVA using ranked variables and analysis of residuals.}
#'   \item{\code{\link{Shirley}}}{Calculates group and interaction effects based on ranks.}
#' }
#'
#' @docType package
#' @name npANCOVA
#' @keywords internal
"_PACKAGE"