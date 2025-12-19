#' Harwell and Serlin Method for Nonparametric ANCOVA
#'
#' Performs the Harwell and Serlin method using ranked response and covariate variables.
#'
#' @param data A data frame containing the variables specified in the formula.
#' @param formula An object of class "formula": a symbolic description of the 
#' model to be fitted. The structure should be `response ~ covariate1 + ... + group`.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{regression_equation}{The summary of the fitted linear model.}
#'   \item{anova}{The ANOVA table from the fitted model.}
#'   \item{statistics}{The Harwell-Serlin test statistic.}
#'   \item{df}{The degrees of freedom for the test.}
#'   \item{p_value}{The p-value of the test.}
#'   \item{data}{The original data frame with added columns for ranks.}
#' }
#'
#' @references
#' Harwell MR, Serlin RCJPB. An empirical study of a proposed test of nonparametric analysis of covariance. 1988;104(2):268.
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
#' # 2. Run the Harwell and Serlin method
#' results <- Harwell_Serlin(
#'   formula = response ~ covariate1 + covariate2 + group,
#'   data = data
#' )
#' 
#' # 3. View the results
#' print(results$p_value)
#' print(paste("Statistic:", results$statistics,"df:", results$df, "P-value:", results$p_value))
#'
#' @importFrom stats lm anova pchisq as.formula terms
#' @export
Harwell_Serlin <- function(data, formula) {

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

  ranked_covariates <- paste0("ranked_", covariates)
  rhs_parts <- c(ranked_covariates, group)
  formula_lm <- as.formula(paste("ranked_response ~", paste(rhs_parts, collapse = " + ")))

  fit <- lm(formula_lm, data = data)
  anova_table <- anova(fit) # Renamed to avoid confusion with the function anova()

  SSB <- anova_table[["Sum Sq"]][length(covariates) + 1]
  SST <- sum(anova_table[["Sum Sq"]])
  df <- length(unique(data[[group]])) - 1

  PSHS <- (nrow(data) - length(covariates) - 1) * (SSB / SST)
  p_value <- 1 - pchisq(PSHS, df)

  return(list(
    regression_equation = summary(fit),
    anova = anova_table,
    statistics = PSHS,
    df = df,
    p_value = p_value,
    data = data
  ))
}