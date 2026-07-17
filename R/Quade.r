#' Quade Method for Nonparametric ANCOVA
#'
#' Performs Quade's ANCOVA using ranked variables and analysis of residuals.
#' The method fits a linear model of the ranked response on the ranked covariates,
#' and then performs an ANOVA on the residuals of that model.
#'
#' @param data A data frame containing the variables specified in the formula.
#' All variables used in the model should be complete; missing values are not allowed.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`.
#' The grouping variable should be the last term in the formula.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{regression_equation}{The summary of the linear model regressing the ranked response on the ranked covariates.}
#'   \item{regression_equation_residuals}{The summary of the model fitting residuals on the group.}
#'   \item{group_means}{A data frame of the mean of residuals for each group.}
#'   \item{group_sds}{A data frame of the standard deviation of residuals for each group.}
#'   \item{anova_summary}{The summary of the ANOVA model performed on the residuals.}
#'   \item{data}{The original data frame augmented with ranked variables and residuals.}
#' }
#'
#' @references
#' Quade, D. (1967). Rank analysis of covariance. \emph{Journal of the American Statistical Association}, 62(320), 1187-1200.
#'
#' Olejnik, S. F., & Algina, J. (1985). A review of nonparametric alternatives to analysis of covariance. \emph{Journal of Educational Research}, 9(1), 51-83.
#'
#' @examples
#' # 1. Create a sample data frame
#' data <- data.frame(
#'   group = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
#'   response = c(16, 60, 82, 126, 137, 44, 67, 87, 100, 142, 17, 28, 105, 149, 160),
#'   covariate1 = c(26, 10, 42, 49, 55, 21, 28, 5, 12, 58, 1, 19, 41, 48, 35),
#'   covariate2 = c(12, 21, 24, 29, 34, 17, 2, 40, 38, 36, 8, 1, 9, 28, 16)
#' )
#'
#' # 2. Run the Quade method
#' results <- Quade(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results)
#' print(results$anova)
#'
#' @importFrom stats lm residuals aov sd as.formula aggregate terms reformulate
#' @export
Quade <- function(formula, data) {

  if (!is.data.frame(data)) {
    stop("The dataset must be a data frame.")
  }
  if (anyNA(data)) {
    stop("The dataset contains missing values.")
  }

  if (grepl("~\\s*\\.$", deparse(formula))) {
    stop("Formula cannot end with '~.'.")
  }

  terms_obj <- terms(formula, data = data)
  response <- all.vars(formula)[1]
  predictors <- attr(terms_obj, "term.labels")

  if (length(response) == 0 || !(response %in% names(data))) {
    stop("No response variable found in the formula or data.")
  }

  group <- predictors[length(predictors)]
  covariates <- predictors[-length(predictors)]

  if (length(covariates) == 0) {
      stop("No covariate found in the formula. At least one is required.")
  }

  if (!is.numeric(data[[response]])) {
    stop("Response must be numeric or integer.")
  }

  if (any(sapply(data[covariates], function(x) !is.numeric(x) && !is.integer(x)))) {
    stop("All covariates must be numeric or integer.")
  }

  if (!is.factor(data[[group]])) {
    data[[group]] <- factor(data[[group]])
  }

  data$ranked_response <- rank(data[[response]], ties.method = "average")

  for (covariate in covariates) {
    data[[paste0("ranked_", covariate)]] <- rank(data[[covariate]])
  }

  ranked_covariates <- paste0("ranked_", covariates)

  formula_lm <- reformulate(termlabels = ranked_covariates, response = "ranked_response")
  lm_fit <- lm(formula_lm, data = data)
  data$residuals <- residuals(lm_fit)

  aov_formula <- reformulate(group, response = "residuals")
  aov_result <- aov(aov_formula, data = data)
  anova_summary <- summary(aov_result)

  lm_fit_residuals <- lm(aov_formula, data = data)

  # Note: `aggregate` can use a formula directly.
  group_mean <- aggregate(aov_formula, data = data, FUN = mean)
  group_sd <- aggregate(aov_formula, data = data, FUN = sd)

  list(
    method_name = "Quade",
    formula = formula,
    regression_equation = summary(lm_fit),
    regression_equation_residuals = summary(lm_fit_residuals),
    group_means = group_mean,
    group_sds = group_sd,
    anova = anova_summary,
    data = data
  )
}