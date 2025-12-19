#' Puri and Sen Method with Unbiased Variance-Covariance Matrix for Nonparametric ANCOVA: One Covariate
#'
#' Performs the Puri and Sen method for a single covariate using an unbiased
#' variance-covariance matrix.
#'
#' @param data A data frame containing the variables specified in the formula.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate + group`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{residuals}{A vector of residuals for each group.}
#'   \item{V}{The unbiased variance-covariance matrix.}
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
#'   covariate1 = c(26, 10, 42, 49, 55, 21, 28, 5, 12, 58, 1, 19, 41, 48, 35)
#' )
#'
#' # 2. Run the Puri and Sen (OU) method
#' results <- Puri_Sen_OU(
#'   formula = response ~ covariate1 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results) 
#' print(paste("Statistic:", results$L_statistic,"df:", results$df, "P-value:", results$p_value))
#'
#' @importFrom stats cor cov pchisq aggregate terms as.formula
#' @export
Puri_Sen_OU <- function(data, formula) {

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

  if (length(predictors) < 2) {
    stop("Formula must include at least one covariate and one group variable.")
  }

  group <- predictors[length(predictors)]
  covariate <- predictors[-length(predictors)]

  if (length(covariate) != 1) {
    stop("The Burnett and Barr method requires exactly one covariate.")
  }

  if (!is.numeric(data[[response]]) && !is.integer(data[[response]])) {
    stop("Response must be numeric or integer.")
  }
  if (!is.numeric(data[[covariate]])) { 
    stop("Covariate must be numeric or integer.")
  }

  if (!is.factor(data[[group]])) {
    data[[group]] <- factor(data[[group]])
  }

  data$ranked_response <- rank(data[[response]])
  data$ranked_covariate <- rank(data[[covariate]])

  f_resp <- as.formula(paste("ranked_response ~", group))
  f_cov <- as.formula(paste("ranked_covariate ~", group))
  response_mean <- aggregate(f_resp, data = data, FUN = mean)
  covariate_mean <- aggregate(f_cov, data = data, FUN = mean)

  r <- cor(data$ranked_response, data$ranked_covariate)
  mean_rank <- (nrow(data) + 1) / 2 # More accurate mean rank

  lvls <- levels(data[[group]])
  response_mean <- response_mean[match(lvls, as.character(response_mean[[group]])), , drop = FALSE]
  covariate_mean <- covariate_mean[match(lvls, as.character(covariate_mean[[group]])), , drop = FALSE]

  k <- length(lvls)

  residuals <- sapply(1:k, function(i) {
    (response_mean$ranked_response[i] - mean_rank) -
      r * (covariate_mean$ranked_covariate[i] - mean_rank)
  })

  V <- cov(data[, c("ranked_response", "ranked_covariate")])
  if (det(V) == 0) {
    warning("Covariance matrix is singular; cannot compute statistic.")
    return(list(L_statistic = NA, p_value = NA))
  }

  inverse_V <- solve(V)

  n_counts_tb <- table(data[[group]])
  n_counts <- as.numeric(n_counts_tb[lvls])

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