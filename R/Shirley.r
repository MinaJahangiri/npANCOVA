#' Shirley Method for Nonparametric ANCOVA
#'
#' Calculates group and interaction effects based on ranked response and
#' covariate variables using changes in R-squared values between models.
#'
#' @param data A data frame containing the variables specified in the formula.
#' All variables used in the model should be complete; missing values are not allowed.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`.
#' The grouping variable should be the last term in the formula.
#'
#' @return A list containing components related to the group and interaction effects, including:
#' \describe{
#'   \item{regression_equation_covariate}{The summary of the model with only covariates.}
#'   \item{regression_equation_covariate_group}{The summary of the model with covariates and group main effects.}
#'   \item{regression_equation_interaction}{The summary of the model including the interaction term.}
#'   \item{statistics_group}{The test statistic for the main group effect.}
#'   \item{df_group}{The degrees of freedom for the group effect.}
#'   \item{p_value_group}{The p-value for the main group effect.}
#'   \item{statistics_interaction}{The test statistic for the interaction effect.}
#'   \item{df_interaction}{The degrees of freedom for the interaction effect.}
#'   \item{p_value_interaction}{The p-value for the interaction effect.}
#'   \item{data}{The original data frame with added rank columns.}
#' }
#'
#' @references
#' Burnett, T. D., & Barr, D. R. (1977). A nonparametric analogy of analysis of covariance. \emph{Educational and Psychological Measurement}, 37(2), 341-348.
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
#' # 2. Run the Shirley method
#' results <- Shirley(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results)
#' print(paste("Statistic:", results$statistics_group,
#'   "df_group:", results$df_group,
#'   "P-value:", results$p_value_group))
#'
#' print(paste("Statistic:", results$statistics_interaction,
#'   "df_interaction:", results$df_interaction,
#'   "P-value:", results$p_value_interaction))
#'
#' @importFrom stats lm pchisq as.formula terms
#' @export
Shirley <- function(formula, data) {

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

  group <- predictors[length(predictors)]
  covariates <- predictors[-length(predictors)]

  if (length(covariates) == 0) {
      stop("No covariate found in the formula. At least one is required.")
  }
  
 if (length(response) == 0 || !(response %in% names(data))) {
    stop("No response variable found in the formula or data.")
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

  data$ranked_response <- rank(data[[response]])

  for (covariate in covariates) {
    data[[paste0("ranked_", covariate)]] <- rank(data[[covariate]])
  }

  ranked_covariates <- paste0("ranked_", covariates)
  covariates_formula <- paste(ranked_covariates, collapse = " + ")

  formula_lm1 <- as.formula(paste("ranked_response ~", covariates_formula))
  formula_lm2 <- as.formula(paste("ranked_response ~", covariates_formula, "+", group))

  interaction_terms <- paste(ranked_covariates, group, sep = ":")
  formula_lm3 <- as.formula(paste(
    "ranked_response ~", covariates_formula, "+", group, "+",
    paste(interaction_terms, collapse = " + ")
  ))

  fit1 <- lm(formula_lm1, data = data)
  fit2 <- lm(formula_lm2, data = data)
  fit3 <- lm(formula_lm3, data = data)

  n <- nrow(data)
  total_mean_square_unadjusted <- n * (n + 1) / 12
  rank_sum_correction <- n * (n^2 - 1) / 12

  r_squared_1 <- summary(fit1)$r.squared
  r_squared_2 <- summary(fit2)$r.squared

  adjusted_sum_squares_group <- (r_squared_2 - r_squared_1) * rank_sum_correction
  statistics_group <- adjusted_sum_squares_group / total_mean_square_unadjusted

  df_group <- length(unique(data[[group]])) - 1
  p_value_group <- 1 - pchisq(statistics_group, df_group)

  r_squared_3 <- summary(fit3)$r.squared
  adjusted_sum_squares_interaction <- (r_squared_3 - r_squared_2) * rank_sum_correction
  statistics_interaction <- adjusted_sum_squares_interaction / total_mean_square_unadjusted

  df_interaction <- df_group * length(covariates)
  p_value_interaction <- 1 - pchisq(statistics_interaction, df_interaction)

  list(
    method_name = "Shirley",
    formula = formula,
    regression_equation_covariate = summary(fit1),
    regression_equation_covariate_group = summary(fit2),
    adjusted_sum_squares_group = adjusted_sum_squares_group,
    total_mean_square_unadjusted = total_mean_square_unadjusted,
    regression_equation_interaction = summary(fit3),
    adjusted_sum_squares_interaction = adjusted_sum_squares_interaction,
    statistics_group = statistics_group,
    df_group = df_group,
    p_value_group = p_value_group,
    statistics_interaction = statistics_interaction,
    df_interaction = df_interaction,
    p_value_interaction = p_value_interaction,
    data = data
  )
}