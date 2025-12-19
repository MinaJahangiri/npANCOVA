#' Burnett and Barr Method for Nonparametric ANCOVA
#'
#' Implements the Burnett and Barr rank-based method for ANCOVA. This method is
#' suitable for models with one response, one covariate, and one grouping variable.
#'
#' @param data A data frame containing the variables specified in the formula.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate + group`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{regression_equation}{The summary of the fitted linear model.}
#'   \item{anova}{The ANOVA table from the fitted model.}
#'   \item{data}{The original data frame with added columns for ranks.}
#' }
#'
#' @references
#' Burnett TD, Barr DRJE, Measurement P. A nonparametric analogy of analysis of covariance. 1977;37(2):341-8.
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
#' # 2. Run the Burnett and Barr method
#' results <- Burnett_Barr(
#'   formula = response ~ covariate1 + group,
#'   data = data
#' )
#'
#' # 3. View the results
#' print(results)
#' print(results$anova)
#'
#' @importFrom stats as.formula lm anova terms
#' @export
Burnett_Barr <- function(data, formula) {
  
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

  if (any(sapply(data[covariate], function(x) !is.numeric(x) && !is.integer(x)))) {
    stop("All covariates must be numeric or integer.")
  }

  if (!is.factor(data[[group]])) {
    data[[group]] <- factor(data[[group]])
  }

  data$ranked_response <- rank(data[[response]])
  data$ranked_covariate <- rank(data[[covariate]])
  data$diff_rank <- data$ranked_response - data$ranked_covariate

  fit_formula <- as.formula(paste("diff_rank ~", group))
  fit <- lm(fit_formula, data = data)

  list(
    regression_equation = summary(fit),
    anova = anova(fit),
    data = data
  )
}