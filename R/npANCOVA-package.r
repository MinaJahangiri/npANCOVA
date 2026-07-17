#' npANCOVA: Nonparametric ANCOVA Methods
#'
#' @description
#' The \pkg{npANCOVA} package provides a unified framework for performing
#' nonparametric Analysis of Covariance (ANCOVA) methods, including Quade,
#' Puri and Sen, McSweeney and Porter, Burnett and Barr, Hettmansperger
#' and McKean, Shirley, and Puri-Sen-Harwell-Serlin. The package provides
#' user-friendly functions to apply these methods in practice.
#'
#' @details
#' This package is designed for researchers and statisticians who need to
#' perform ANCOVA when the assumptions of classical parametric ANCOVA,
#' such as normality of residuals is not met, or when the response variable is ordinal.
#'
#' @section Primary Functions:
#' \describe{
#'   \item{\code{\link{Burnett_Barr}}}{The Burnett and Barr rank-based method compares groups while adjusting for only one covariate.}
#'   \item{\code{\link{Harwell_Serlin}}}{Performs the Harwell and Serlin method using ranked response and covariate variables.}
#'   \item{\code{\link{Hettmansperger_McKean}}}{Applies rank-based residual analysis for ANCOVA, and then performs an ANOVA
#'    on the (weighted) ranked residuals.}
#'   \item{\code{\link{McSweeney_Porter}}}{Performs rank-based ANCOVA with and without an interaction term between the covariates and 
#'    the group.}
#'   \item{\code{\link{Puri_Sen_UVCM}}}{Performs the Puri and Sen method using an unbiased variance-covariance matrix.}
#'   \item{\code{\link{Puri_Sen_BVCM}}}{Performs the Puri and Sen method using a biased variance-covariance matrix.}
#'   \item{\code{\link{Quade}}}{Performs Quade’s ANCOVA using ranked variables and analysis of residuals using ANOVA.}
#'   \item{\code{\link{Shirley}}}{Calculates group and interaction effects based on ranked response and covariate variables using
#'    changes in R-squared values between models.}
#'   \item{\code{\link{summary_npANCOVA}}}{Formats and prints the inferential results obtained from a fitted nonparametric ANCOVA method.}
#'   \item{\code{\link{pairwise_npANCOVA}}}{Performs pairwise post hoc comparisons based on the chosen nonparametric ANCOVA method.}
#' }
#'
#' @docType package
#' @name npANCOVA-package
#' @aliases npANCOVA-package
#' @keywords internal
#' @importFrom stats aggregate anova aov as.formula cor cov cov.wt lm p.adjust pchisq qchisq reformulate residuals sd terms
#' @importFrom utils combn head
"_PACKAGE"