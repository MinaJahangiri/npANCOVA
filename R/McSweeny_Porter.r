#' McSweeny and Porter Method for Nonparametric ANCOVA
#'
#' Performs rank-based ANCOVA with and without an interaction term between the covariates and the group.
#'
#' @param data A data frame containing the variables specified in the formula.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{regression_equation_covariate}{Summary of the model with only covariates.}
#'   \item{regression_equation_covariate_group}{Summary of the model with covariates and group main effects.}
#'   \item{group_effect}{The result of an ANOVA test for group effect.}
#'   \item{interaction_effect}{The result of an ANOVA test for interaction effect between group and covariate variables.}
#'   \item{regression_equation_interaction}{Summary of the model including the interaction term.}
#'   \item{data}{The original data frame with added columns for ranks.}
#' }
#'
#' @references
#' McSweeney M, Porter AJOp. Small sample properties of nonparametric index of response and rank analysis of covariance. 1971;16.
#'
#' Olejnik SF, Algina JJER. A review of nonparametric alternatives to analysis of covariance. 1985;9(1):51-83.
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
#' # 2. Run the McSweeny and Porter method
#' results <- McSweeny_Porter(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results)
#' print(results$group_effect)
#' print(results$interaction_effect)
#'
#' @importFrom stats lm anova as.formula terms
#' @export
McSweeny_Porter <- function(data, formula) {
  if (!is.data.frame(data)) {
    stop("Error: The dataset must be a data frame.")
  }
  if (anyNA(data)) {
    stop("Error: The dataset contains missing values.")
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
    stop("Error: The response variable must be numeric or integer.")
  }

  if (any(sapply(data[covariates], function(x) !is.numeric(x) && !is.integer(x)))) {
    stop("Error: All covariates must be numeric or integer.")
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
  formula_lm3 <- as.formula(paste("ranked_response ~",
    paste(c(covariates_formula, group, paste(interaction_terms, collapse = " + ")),
      collapse = " + "
    )
  ))
  fit1 <- lm(formula_lm1, data = data)
  fit2 <- lm(formula_lm2, data = data)
  group_effect <- anova(fit1, fit2)

  fit3 <- lm(formula_lm3, data = data)
  interaction_effect <- anova(fit2, fit3)

  list(
    regression_equation_covariate = summary(fit1),
    regression_equation_covariate_factor = summary(fit2),
    group_effect = group_effect,
    interaction_effect = interaction_effect,
    regression_equation_interaction = summary(fit3),
    data = data
  )
}