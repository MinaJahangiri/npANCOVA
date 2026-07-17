npANCOVA <- function(formula, data,
              method = c("Burnett_Barr",
                         "Harwell_Serlin",
                         "Puri_Sen_UVCM",
                         "Puri_Sen_BVCM",
                         "Quade",
                         "Shirley",
                         "McSweeney_Porter",
                         "Hettmansperger_McKean")) {

  method <- match.arg(method)

  methods_list <- list(
    Burnett_Barr        = Burnett_Barr,
    Harwell_Serlin        = Harwell_Serlin,
    Puri_Sen_UVCM           = Puri_Sen_UVCM,
    Puri_Sen_BVCM           = Puri_Sen_BVCM,
    Quade                 = Quade,
    Shirley               = Shirley,
    McSweeney_Porter       = McSweeney_Porter,
    Hettmansperger_McKean = Hettmansperger_McKean
  )

  selected_fun <- methods_list[[method]]
  res <- selected_fun(data = data, formula = formula)
  res$method_name <- method
  res$formula <- formula
  class(res) <- "npANCOVA"
  return(res)
}

##############################################

extract_npancova_pvalue <- function(fit) {
  method <- fit$method_name
  
  if (is.null(method)) {
    stop("method_name not found in fitted object.")
  }
  
  pval <- switch(
    method,
    
    "Burnett_Barr" = {
      a <- fit$anova
      if (!is.null(a) && "Pr(>F)" %in% colnames(a)) {
        a[1, "Pr(>F)"]
      } else {
        NA_real_
      }
    },
    
    "Harwell_Serlin" = {
      fit$p_value
    },
    
    "Puri_Sen_UVCM" = {
      fit$p_value
    },
    
    "Puri_Sen_BVCM" = {
      fit$p_value
    },
    
    "Quade" = {
      a <- fit$anova[[1]]
      if (!is.null(a) && "Pr(>F)" %in% colnames(a)) {
        a[1, "Pr(>F)"]
      } else {
        NA_real_
      }
    },
    
    "Shirley" = {
      fit$p_value_group
    },
    
    "McSweeney_Porter" = {
      a <- fit$group_effect
      if (!is.null(a) && "Pr(>F)" %in% colnames(a)) {
        a[2, "Pr(>F)"]
      } else {
        NA_real_
      }
    },
    
    "Hettmansperger_McKean" = {
      if (!is.null(fit) && "p_value" %in% names(fit)) {
        as.numeric(fit$p_value)
      } else {
        NA_real_
      }
    },
    stop(paste("Unsupported method:", method))
  )
  
  as.numeric(pval)
}
#' Pairwise Post Hoc Comparisons for Nonparametric ANCOVA
#'
#' @description
#' Performs pairwise post hoc comparisons for nonparametric ANCOVA
#' methods by fitting the selected nonparametric method separately to each
#' pair of groups. For every pairwise comparison, the function extracts the raw
#' p-value from the fitted model and computes multiple-testing adjusted p-values.
#'
#' @details
#' All possible pairwise combinations of group levels are generated, and nonparametric
#' ANCOVA method is fitted to each two-group subset using the selected method.
#'
#' Raw p-values are extracted from each method and adjusted using one or more
#' methods supported by \code{\link[stats]{p.adjust}}. The adjusted p-values are
#' returned together with significance decisions based on the specified significance
#' level \code{alpha}.
#'
#' Supported methods are:
#' \code{"Burnett_Barr"},
#' \code{"Harwell_Serlin"},
#' \code{"Puri_Sen_UVCM"},
#' \code{"Puri_Sen_BVCM"},
#' \code{"Quade"},
#' \code{"Shirley"},
#' \code{"McSweeney_Porter"},
#' and \code{"Hettmansperger_McKean"}.
#'
#' @param formula An object of class "formula": a symbolic description of the model to be fitted.
#' The structure should be ‘response ~ covariate1 + ... + group‘. 
#' The grouping variable should be the last term in the formula.
#'
#' @param data A data frame containing the variables specified  in the model. All variables 
#' used in the model should be complete; missing values are not allowed.
#' @param method A character string specifying the nonparametric ANCOVA method.
#' Should be one of \code{"Burnett_Barr"}, \code{"Harwell_Serlin"},
#' \code{"Puri_Sen_UVCM"}, \code{"Puri_Sen_BVCM"}, \code{"Quade"}, \code{"Shirley"},
#' \code{"McSweeney_Porter"}, or \code{"Hettmansperger_McKean"}.
#' @param p_adjust_methods A character vector specifying one or more p-value
#' adjustment methods to be applied to the raw pairwise p-values using
#' \code{\link[stats]{p.adjust}}.
#' @param alpha A numeric significance level used to determine whether adjusted
#' p-values are statistically significant.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{method}{The selected nonparametric ANCOVA method.}
#'   \item{formula}{The model formula used in the analysis.}
#'   \item{alpha}{The significance level used for adjusted comparisons.}
#'   \item{comparisons}{A data frame containing pairwise group comparisons,
#'   raw (unadjusted) p-values, adjusted p-values, and significance decisions.}
#' }
#'
#' @examples
#' data <- data.frame(
#'   group = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
#'   response = c(16, 60, 82, 126, 137, 44, 67, 87, 100, 142, 17, 28, 105, 149, 160),
#'   covariate1 = c(26, 10, 42, 49, 55, 21, 28, 5, 12, 58, 1, 19, 41, 48, 35)
#' )
#'
#' # Burnett_Barr
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Burnett_Barr"
#' )
#'
#' # Harwell_Serlin
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Harwell_Serlin"
#' )
#'
#' # Puri_Sen_UVCM
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Puri_Sen_UVCM"
#' )
#'
#' # Puri_Sen_BVCM
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Puri_Sen_BVCM"
#' )
#'
#'
#' # Quade
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Quade"
#' )
#'
#' # Shirley
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Shirley"
#' )
#'
#' # McSweeney_Porter
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "McSweeney_Porter"
#' )
#'
#' # Hettmansperger_McKean
#' pairwise_npANCOVA(
#'   response ~ covariate1 + group,
#'   data = data,
#'   method = "Hettmansperger_McKean"
#' )
#'
#' @export
#'
pairwise_npANCOVA <- function(formula, data,
                              method = c("Burnett_Barr",
                                         "Harwell_Serlin",
                                         "Puri_Sen_UVCM",
                                         "Puri_Sen_BVCM",
                                         "Quade",
                                         "Shirley",
                                         "McSweeney_Porter",
                                         "Hettmansperger_McKean"),
                              p_adjust_methods = c("holm", "hochberg", "hommel",
                                                   "bonferroni", "BH", "BY", "fdr"),
                              alpha = 0.05) {
  
  method <- match.arg(method)
  
  if (!is.data.frame(data)) {
    stop("The dataset must be a data frame.")
  }
  
  if (grepl("~\\s*\\.$", deparse(formula))) {
    stop("Formula cannot end with '~.'.")
  }
  
  terms_obj <- terms(formula, data = data)
  response <- all.vars(formula)[1]
  predictors <- attr(terms_obj, "term.labels")
  
  if (length(response) == 0 || !(response %in% names(data))) {
    stop("No response variable found in formula or data.")
  }
  
  if (length(predictors) < 2) {
    stop("Formula must contain at least one covariate and one group variable.")
  }
  
  group_var <- predictors[length(predictors)]
  
  if (!(group_var %in% names(data))) {
    stop("Group variable not found in data.")
  }
  
  if (!is.factor(data[[group_var]])) {
    data[[group_var]] <- factor(data[[group_var]])
  }
  
  group_levels <- levels(droplevels(data[[group_var]]))
  
  if (length(group_levels) < 2) {
    stop("At least two groups are required.")
  }
  
  comps <- combn(group_levels, 2, simplify = FALSE)
  
  raw_pvals <- numeric(length(comps))
  comparison_labels <- character(length(comps))
  fit_list <- vector("list", length(comps))
  
  for (i in seq_along(comps)) {
    g1 <- comps[[i]][1]
    g2 <- comps[[i]][2]
    
    sub_data <- data[data[[group_var]] %in% c(g1, g2), , drop = FALSE]
    sub_data[[group_var]] <- droplevels(factor(sub_data[[group_var]]))
    
    fit_i <- npANCOVA(
      data = sub_data,
      formula = formula,
      method = method
    )
    
    p_i <- extract_npancova_pvalue(fit_i)
    
    comparison_labels[i] <- paste(g1, "vs", g2)
    raw_pvals[i] <- p_i
    fit_list[[i]] <- fit_i
  }
  
  out <- data.frame(
    Comparison = comparison_labels,
    Raw_P_Value = raw_pvals,
    stringsAsFactors = FALSE
  )
  
  for (padj_method in p_adjust_methods) {
    adj_p <- p.adjust(raw_pvals, method = padj_method)
    
    method_label <- toupper(padj_method)
    if (padj_method == "bonferroni") method_label <- "Bonferroni"
    if (padj_method == "holm") method_label <- "Holm"
    if (padj_method == "hochberg") method_label <- "Hochberg"
    if (padj_method == "hommel") method_label <- "Hommel"
    if (padj_method == "BH") method_label <- "BH"
    if (padj_method == "BY") method_label <- "BY"
    if (padj_method == "fdr") method_label <- "FDR"
    
    out[[paste0("P_adj_", method_label)]] <- adj_p
    out[[paste0("Significant_", method_label)]] <- ifelse(adj_p < alpha, "Yes", "No")
  }
  
  result <- list(
    method = method,
    formula = formula,
    alpha = alpha,
    comparisons = out
  )
  
  return(result)
}
