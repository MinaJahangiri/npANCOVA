#' Puri and Sen Method with Biased Variance-Covariance Matrix for Nonparametric ANCOVA: Multiple Covariates
#'
#' Performs the Puri and Sen method for multiple covariates using a biased
#' variance-covariance matrix.
#'
#' @param data A data frame containing the variables specified in the formula.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{residuals}{A vector of residuals for each group.}
#'   \item{V}{The biased variance-covariance matrix.}
#'   \item{inverse_V}{The inverse of the variance-covariance matrix.}
#'   \item{L_statistic}{The Puri and Sen L-statistic.}
#'   \item{df}{The degrees of freedom for the test.}
#'   \item{p_value}{The corresponding p-value of the L-statistic.}
#'   \item{data}{The original data frame with added columns for ranks.}
#' }
#'
#' @references
#' Puri ML, Sen PKJAoMS. Analysis of covariance based on general rank scores. 1969;40(2):610-8.
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
#' # 2. Run the Puri and Sen (MB) method
#' results <- Puri_Sen_MB(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results) 
#' print(paste("Statistic:", results$L_statistic,"df:", results$df, "P-value:", results$p_value))
#'
#' @importFrom stats lm cov.wt pchisq terms as.formula
#' @export
Puri_Sen_MB <- function(data, formula) {
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

  data$ranked_response <- rank(data[[response]])

  for (covariate in covariates) {
    data[[paste0("ranked_", covariate)]] <- rank(data[[covariate]])
  }

  unique_factors <- levels(data[[group]])
  factor_mean <- data.frame(factor = unique_factors)

  for (covariate in covariates) {
    factor_mean[[paste0("ranked_", covariate, "_mean")]] <- sapply(unique_factors, function(f) {
      mean(data[data[[group]] == f, paste0("ranked_", covariate)], na.rm = TRUE)
    })
  }

  factor_mean$ranked_response_mean <- sapply(unique_factors, function(f) {
    mean(data[data[[group]] == f, "ranked_response"], na.rm = TRUE)
  })

  ranked_cols <- c("ranked_response", paste0("ranked_", covariates))
  mean_rank <- mean(rowMeans(data[ranked_cols])) # Simplified mean rank calculation

  factor_mean_response <- factor_mean$ranked_response_mean - mean_rank

  factor_mean_covariate <- matrix(NA, nrow = length(covariates), ncol = length(unique_factors))
  for (i in seq_along(unique_factors)) {
    factor_data <- data[data[[group]] == unique_factors[i], paste0("ranked_", covariates), drop = FALSE]
    factor_mean_covariate[, i] <- colMeans(factor_data, na.rm = TRUE) - mean_rank
  }

  formula_lm <- as.formula(paste("ranked_response ~", paste(paste0("ranked_", covariates), collapse = "+")))
  lm_fit <- lm(formula_lm, data = data)

  residuals <- numeric(length(unique_factors))
  coef_values <- summary(lm_fit)$coef[-1, 1]

  for (i in seq_along(unique_factors)) {
    residuals[i] <- factor_mean_response[i] - sum(coef_values * factor_mean_covariate[, i])
  }

  V <- cov.wt(data[ranked_cols], method = "ML")$cov
  if (det(V) == 0) {
    warning("Covariance matrix is singular; cannot compute statistic.")
    return(list(L_statistic = NA, p_value = NA))
  }
  inverse_V <- solve(V)

  k <- length(unique_factors)
  n_counts_tb <- table(data[[group]])
  n_counts <- as.numeric(n_counts_tb[unique_factors])

  L <- sum(n_counts * residuals^2)
  L_statistics <- L * inverse_V[1, 1]
  df <- k - 1
  p_value <- 1 - pchisq(L_statistics, df)

  return(list(
    residuals = residuals,
    V = V,
    inverse_V = inverse_V,
    L_statistic = L_statistics,
    df = df,
    p_value = p_value,
    data = data
  ))
}