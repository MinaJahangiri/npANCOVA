#' Summarize Nonparametric ANCOVA Results
#'
#' This function formats and prints inferential results obtained from a fitted
#' nonparametric ANCOVA object. 
#'
#' @param object A fitted result object returned by one
#'   of the supported methods: \code{Burnett_Barr()},
#'   \code{Harwell_Serlin()}, \code{Hettmansperger_McKean()},
#'   \code{McSweeney_Porter()}, \code{Puri_Sen_BVCM()},
#'   \code{Puri_Sen_UVCM()}, \code{Quade()}, or \code{Shirley()}.
#'
#' @details
#' The function uses the \code{method_name} stored in \code{object} and prints the
#' appropriate inferential summary according to the detected nonparametric
#' ANCOVA method.
#'
#' \itemize{
#'   \item For \code{"Burnett_Barr"} and \code{"Quade"},
#'   the function prints the F-statistic, numerator and denominator degrees of freedom, and p-value.
#'
#'   \item For \code{"Hettmansperger_McKean"}, \code{"Harwell_Serlin"}, \code{"Puri_Sen_BVCM"}, and
#'   \code{"Puri_Sen_UVCM"}, the function prints the chi-squared test statistic, degrees
#'   of freedom, and p-value.
#'
#'   \item For \code{"McSweeney_Porter"} and \code{"Shirley"}, both group and
#'   interaction effects are printed when available.
#' }
#'
#' @return Prints a formatted summary of the fitted object containing the following components:
#' \describe{
#'   \item{statistics}{The test statistic.}
#'   \item{df}{The degrees of freedom for the test.}
#'   \item{p_value}{The p-value of the test.}
#' } 
#'
#' @seealso
#' \code{\link{Burnett_Barr}},
#' \code{\link{Quade}},
#' \code{\link{Hettmansperger_McKean}},
#' \code{\link{Harwell_Serlin}},
#' \code{\link{McSweeney_Porter}},
#' \code{\link{Shirley}},
#' \code{\link{Puri_Sen_BVCM}},
#' \code{\link{Puri_Sen_UVCM}},
#'
#' @examples
#' data <- data.frame(
#'   group = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
#'   response = c(16, 60, 82, 126, 137, 44, 67, 87, 100, 142, 17, 28, 105, 149, 160),
#'   covariate1 = c(26, 10, 42, 49, 55, 21, 28, 5, 12, 58, 1, 19, 41, 48, 35),
#'   covariate2 = c(12, 21, 24, 29, 34, 17, 2, 40, 38, 36, 8, 1, 9, 28, 16)
#' )
#'
#' formula_one_cov <- response ~ covariate1 + group
#' formula_two_cov <- response ~ covariate1 + covariate2 + group
#'
#' # Method requiring exactly one covariate
#' fit_bb <- Burnett_Barr(formula_one_cov, data = data)
#' summary_npANCOVA(fit_bb)
#'
#' # Methods allowing multiple covariates
#' fit_hs <- Harwell_Serlin(formula_two_cov, data = data)
#' summary_npANCOVA(fit_hs)
#'
#' fit_hm <- Hettmansperger_McKean(formula_two_cov, data = data)
#' summary_npANCOVA(fit_hm)
#'
#' fit_mp <- McSweeney_Porter(formula_two_cov, data = data)
#' summary_npANCOVA(fit_mp)
#'
#' fit_psmb <- Puri_Sen_BVCM(formula_two_cov, data = data)
#' summary_npANCOVA(fit_psmb)
#'
#' fit_psmu <- Puri_Sen_UVCM(formula_two_cov, data = data)
#' summary_npANCOVA(fit_psmu)
#'
#' fit_quade <- Quade(formula_two_cov, data = data)
#' summary_npANCOVA(fit_quade)
#'
#' fit_shirley <- Shirley(formula_two_cov, data = data)
#' summary_npANCOVA(fit_shirley)
#'
#' @export
summary_npANCOVA <- function(object) {

  detect_method <- function(x) {
    method <- x$method_name

    valid_methods <- c(
      "Burnett_Barr",
      "Harwell_Serlin",
      "Puri_Sen_UVCM",
      "Puri_Sen_BVCM",
      "Quade",
      "Shirley",
      "McSweeney_Porter",
      "Hettmansperger_McKean"
    )

    if (is.null(method) || !is.character(method) || length(method) != 1L || !nzchar(method)) {
      stop(
        "method_name not found in fitted object. ",
        "Please fit the model using npANCOVA() or store method_name in the returned object."
      )
    }

    if (!(method %in% valid_methods)) {
      stop("Unsupported method_name in fitted object: ", method)
    }

    method
  }

  get_anova_values <- function(anova_obj) {
    if (is.list(anova_obj) && !is.data.frame(anova_obj)) {
      anova_obj <- anova_obj[[1]]
    }

    f_val <- anova_obj[1, "F value"]
    df_num <- anova_obj[1, "Df"]

    row_names <- trimws(rownames(anova_obj))
    res_idx <- which(row_names == "Residuals")

    if (length(res_idx) > 0) {
      df_den <- anova_obj[res_idx[1], "Df"]
    } else {
      df_den <- anova_obj[nrow(anova_obj), "Df"]
    }

    p_val <- anova_obj[1, "Pr(>F)"]

    list(F = f_val, Df_num = df_num, Df_den = df_den, P = p_val)
  }

  get_signif_stars <- function(p) {
    if (length(p) == 0 || is.null(p) || is.na(p)) {
      return("")
    }
    if (p < 0.001) {
      return("***")
    }
    if (p < 0.01) {
      return("**")
    }
    if (p < 0.05) {
      return("*")
    }
    if (p < 0.1) {
      return(".")
    }
    ""
  }

  format_numeric <- function(x, digits = 4) {
    if (length(x) == 0 || is.null(x) || is.na(x)) {
      return("NA")
    }
    if (is.numeric(x)) {
      return(formatC(x, digits = digits, format = "f"))
    }
    as.character(x)
  }

  format_p_value <- function(p, digits = 4) {
    if (length(p) == 0 || is.null(p) || is.na(p)) {
      return("NA")
    }

    out <- if (p < 0.0001) {
      "< 0.0001"
    } else {
      formatC(p, digits = digits, format = "f")
    }

    stars <- get_signif_stars(p)

    if (nzchar(stars)) {
      paste0(out, " ", stars)
    } else {
      out
    }
  }

  get_formula_text <- function(x) {
    if (!is.null(x$formula)) {
      return(paste(deparse(x$formula), collapse = " "))
    }

    candidates <- c(
      "regression_equation",
      "regression_equation_covariate",
      "regression_equation_residuals"
    )

    for (cand in candidates) {
      if (!is.null(x[[cand]]) && !is.null(x[[cand]]$call$formula)) {
        return(paste(deparse(x[[cand]]$call$formula), collapse = " "))
      }
    }

    "Not specified"
  }

  line_width <- 58
  rule_top <- paste(rep("=", line_width), collapse = "")
  rule_mid <- paste(rep("-", line_width), collapse = "")
  label_width <- 20

  print_field <- function(label, value) {
    cat(sprintf("  %-*s : %s\n", label_width, label, value))
  }

  print_section <- function(title) {
    cat(sprintf("\n[%s]\n", title))
  }

  method_name <- detect_method(object)
  formula_text <- get_formula_text(object)

  cat("\n", rule_top, "\n", sep = "")
  cat(" Nonparametric ANCOVA (npANCOVA)\n")
  cat(rule_top, "\n", sep = "")
  print_field("Method", method_name)
  print_field("Formula", formula_text)
  cat(rule_mid, "\n", sep = "")

  has_stars <- FALSE

  switch(
    method_name,

    "Burnett_Barr" = {
      vals <- get_anova_values(object$anova)
      print_section("Group Effect")
      print_field("Statistic", format_numeric(vals$F))
      print_field("Degrees of Freedom", sprintf("%s, %s", vals$Df_num, vals$Df_den))
      print_field("P-value", format_p_value(vals$P))
      has_stars <- nzchar(get_signif_stars(vals$P))
    },

    "Quade" = {
      vals <- get_anova_values(object$anova)
      print_section("Group Effect")
      print_field("Statistic", format_numeric(vals$F))
      print_field("Degrees of Freedom", sprintf("%s, %s", vals$Df_num, vals$Df_den))
      print_field("P-value", format_p_value(vals$P))
      has_stars <- nzchar(get_signif_stars(vals$P))
    },

    "Hettmansperger_McKean" = {
      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$statistics))
      print_field("Degrees of Freedom", object$df)
      print_field("P-value", format_p_value(object$p_value))
      has_stars <- nzchar(get_signif_stars(object$p_value))
    },

    "Harwell_Serlin" = {
      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$statistics))
      print_field("Degrees of Freedom", object$df)
      print_field("P-value", format_p_value(object$p_value))
      has_stars <- nzchar(get_signif_stars(object$p_value))
    },

    "McSweeney_Porter" = {
      group_p <- object$group_effect[2, "Pr(>F)"]
      interaction_p <- object$interaction_effect[2, "Pr(>F)"]

      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$group_effect[2, "F"]))
      print_field(
        "Degrees of Freedom",
        sprintf("%s, %s", object$group_effect[2, "Df"], object$group_effect[2, "Res.Df"])
      )
      print_field("P-value", format_p_value(group_p))

      print_section("Interaction Effect")
      print_field("Statistic", format_numeric(object$interaction_effect[2, "F"]))
      print_field(
        "Degrees of Freedom",
        sprintf("%s, %s", object$interaction_effect[2, "Df"], object$interaction_effect[2, "Res.Df"])
      )
      print_field("P-value", format_p_value(interaction_p))

      has_stars <- nzchar(get_signif_stars(group_p)) ||
        nzchar(get_signif_stars(interaction_p))
    },

    "Puri_Sen_UVCM" = {
      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$L_statistic))
      print_field("Degrees of Freedom", object$df)
      print_field("P-value", format_p_value(object$p_value))
      has_stars <- nzchar(get_signif_stars(object$p_value))
    },

    "Puri_Sen_BVCM" = {
      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$L_statistic))
      print_field("Degrees of Freedom", object$df)
      print_field("P-value", format_p_value(object$p_value))
      has_stars <- nzchar(get_signif_stars(object$p_value))
    },

    "Shirley" = {
      print_section("Group Effect")
      print_field("Statistic", format_numeric(object$statistics_group))
      print_field("Degrees of Freedom", object$df_group)
      print_field("P-value", format_p_value(object$p_value_group))

      print_section("Interaction Effect")
      print_field("Statistic", format_numeric(object$statistics_interaction))
      print_field("Degrees of Freedom", object$df_interaction)
      print_field("P-value", format_p_value(object$p_value_interaction))

      has_stars <- nzchar(get_signif_stars(object$p_value_group)) ||
        nzchar(get_signif_stars(object$p_value_interaction))
    }
  )

  if (has_stars) {
    cat("\n", rule_mid, "\n", sep = "")
    cat(" Signif. codes:  0 '***'  0.001 '**'  0.01 '*'  0.05 '.'  0.1 ' '  1\n", sep = "")
  }

  cat(rule_top, "\n", sep = "")
}
