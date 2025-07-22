#' Tidy the Output of Model Fit
#'
#' Takes the output of \code{\link{Fit_piModel}} and returns it into a user-friendly formats.
#'
#' @param fit A list returned by \code{\link{Fit_piModel}}, which includes model parameters,
#' statistical metrics, and (optionally) a Hessian matrix.
#' @param clevel Numeric. Confidence level for the intervals (e.g., 0.95).
#' If `NULL`, confidence intervals are skipped unless available from the fit.
#' Defaults to `0.95` if Hessian is present and user does not specify it.
#'
#' @return A named list with the following elements:
#' \itemize{
#'
#'   \item{\code{model}: }{The model name or expression used for fitting.}
#'   \item{\code{params}: }{A tibble of estimated parameters, standard errors, and confidence intervals (if available).}
#'   \item{\code{stats_metric}: }{A tibble of fit quality metrics, including SQA and information criteria.}
#' }
#' @examples
#' df <- tibble::tibble(
#'   I = c(0, 33, 67, 100, 133, 167, 200),
#'   P = c(0, 15, 19, 20, 20, 20, 20)
#' )
#'
#' fit <- Fit_piModel(df)
#' tidy_fit <- Tidy_piCurve(fit)
#'
#'
#' @export
Tidy_piCurve <- function(fit, clevel = NULL) {
    has_hessian <- !is.null(fit$hessian)

    # Default confidence level handling
    if (has_hessian && is.null(clevel)) {
        clevel <- 0.95
    }

    # Parameter table initialization
    tb_par <- tibble::tibble(params = names(fit$par), Est_Val = fit$par)

    # Add confidence intervals if available
    if (has_hessian && !is.null(clevel)) {
        fit_ci <- ConfInt_piCurve(fit, Clevel = clevel)

        # Generate dynamic column names
        clevel_suffix <- as.character(clevel*100)
        lb_col_name <- paste0("LB_", clevel_suffix, "%")
        ub_col_name <- paste0("UB_", clevel_suffix, "%")

        tb_par$Est_Err <- fit_ci[startsWith(names(fit_ci), "Est_Err")]
        tb_par[[lb_col_name]] <- fit_ci[startsWith(names(fit_ci), "L_CI")]
        tb_par[[ub_col_name]] <- fit_ci[startsWith(names(fit_ci), "U_CI")]
    }

    # Stats table
    stat_names <- names(fit$SQA)
    stat_vals <- fit$SQA
    if (!is.null(fit$info_criteria)) {
        stat_names <- c(stat_names, names(fit$info_criteria))
        stat_vals <- c(stat_vals, fit$info_criteria)
    }
    tb_stat <- tibble::tibble(metrics = stat_names, values = stat_vals)

    # Output
    out <- list(
        model = fit$model,
        params = tb_par,
        stats_metric = tb_stat
    )

    return(out)
}
