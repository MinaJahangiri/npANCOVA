#' Hettmansperger and McKean Method for ANCOVA
#'
#' Applies rank-based residual analysis for ANCOVA. This method involves
#' fitting a model of the response on the covariate (s), calculating residuals,
#' ranking them, and then performing an ANOVA on the (weighted) ranked residuals.
#'
#' @param data A data frame containing the variables specified in the formula.
#' All variables used in the model should be complete; missing values are not allowed.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`. 
#' The grouping variable should be the last term in the formula.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{regression_equation_covariate}{The summary of the initial model fitting response on covariates.}
#'   \item{regression_equation_residuals}{The summary of the model fitting weighted ranked residuals on the group.}
#'   \item{group_means}{A data frame of the mean of weighted ranked residuals for each group.}
#'   \item{group_sds}{A data frame of the standard deviation of weighted ranked residuals for each group.}
#'   \item{anova}{The ANOVA table for the model based on weighted ranked residuals.}
#'   \item{statistic}{The test statistic.}
#'   \item{df}{The degrees of freedom for the test.}
#'   \item{p_value}{The p-value of the test.}
#'   \item{data}{The original data frame augmented with residuals, ranked residuals, and weighted ranked residuals.}
#' }
#'
#' @references
#' Hettmansperger, T. P., & McKean, J. W. (1977). A robust alternative based on ranks to least squares in analyzing linear models. \emph{Technometrics}, 19(3), 275-284.
#'
#' Hettmansperger, T. P., & McKean, J. W. (1983). A geometric interpretation of inferences based on ranks in the linear model. \emph{Journal of the American Statistical Association}, 78(384), 885-893.
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
#' # 2. Run the Hettmansperger and McKean method
#' results <- Hettmansperger_McKean(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results)
#' print(paste("Statistic:", results$statistics,
#' "df:", results$df,
#' "P-value:", results$p_value))
#'
#' @importFrom stats lm residuals anova aggregate sd terms reformulate
#' @export
Hettmansperger_McKean <- function(formula, data) {

  if (!is.data.frame(data)) {
    stop("The dataset must be a data frame.")
  }
  if (any(is.na(data))) {
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

  if (!is.numeric(data[[response]]) && !is.integer(data[[response]])) {
    stop("Response must be numeric or integer.")
  }

  if (any(sapply(data[covariates], function(x) !is.numeric(x) && !is.integer(x)))) {
    stop("All covariates must be numeric or integer.")
  }

  if (!is.factor(data[[group]])) {
    data[[group]] <- factor(data[[group]])
  }

  formula_lm1 <- reformulate(termlabels = covariates, response = response)
  fit <- lm(formula_lm1, data = data)

  data$residuals <- residuals(fit)
  data$rank_residuals <- rank(data$residuals, ties.method = "average")
  data$weighted_rank_residuals <- sqrt(12) *
    ((data$rank_residuals / (nrow(data) + 1)) - 0.5)

  formula_lm2 <- reformulate(termlabels = group, response = "weighted_rank_residuals")
  fit_residuals <- lm(formula_lm2, data = data)
  anova_table <- anova(fit_residuals)

  group_row <- anova_table[1, , drop = FALSE]

  statistics <- unname(group_row[["Sum Sq"]])
  df_group <- unname(group_row[["Df"]])
  p_value <- pchisq(statistics, df = df_group, lower.tail = FALSE)
  critical_value <- qchisq(0.95, df = df_group)

  group_mean <- aggregate(formula_lm2, data = data, FUN = mean)
  group_sd <- aggregate(formula_lm2, data = data, FUN = sd)

  result <- list(
    method_name = "Hettmansperger_McKean",
    formula = formula,
    regression_equation_covariate = summary(fit),
    regression_equation_residuals = summary(fit_residuals),
    group_means = group_mean,
    group_sds = group_sd,
    anova = anova_table,
    statistics = statistics,
    df = df_group,
    p_value = p_value,
    data = data
  )
}


