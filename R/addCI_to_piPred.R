#' Add Confidence Interval to PI Predictions
#'
#' Generates confidence intervals for predicted photosynthesis values (`P`) given a fitted PI curve model.
#'
#' @param Fitted_Model A model object returned by \code{\link{Fit_piModel}} including parameter estimates and model name.
#' @param irrad A numeric vector of irradiance (`I`) values at which to compute prediction intervals.
#' @param n_cores Integer. Number of CPU cores to use for parallel sampling (default: 2).
#' @param n_sample Integer. Number of parameter samples to draw from the confidence intervals (default: 1000).
#' @param Clevel Confidence level for interval computation (default: 0.95).
#'
#' @return A data frame with columns:
#' \itemize{
#'   \item `irradiance`: input irradiance values
#'   \item `ymin`: lower bound of predicted `P` at each `I`
#'   \item `ymax`: upper bound of predicted `P` at each `I`
#' }
#'
#' @details
#' This function first extracts confidence intervals for each estimated model parameter using \code{\link{ConfInt_piCurve}}.
#' It then samples `n_sample` parameter sets from these intervals and generates predictions for each.
#' Finally, it computes the confidence interval bounds of predicted photosynthetic rate at each irradiance level.
#'
#' Parallel computation is used to speed up sampling. The interval is two-tailed and centered based on the requested confidence level.
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows group_by summarise
#' @importFrom parallel mclapply
#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#' df$P <- Model_piCurve(parameters = params, model_name = "LS5", data = df) +
#'   5 * rnorm(25, 0, 0.25)
#'
#' fit <- Fit_piModel(parameters = c(params, StDev = 2), STATapp = "MLE",
#'                    Hessian = TRUE, data = df)
#'
#' cal_ci <- addCI_to_piPred(fit, irrad = df$I, n_cores = 2, n_sample = 500)
#'
#' df_ci <- tibble::tibble(
#'   cal_ci,
#'   .obs = df$P,
#'   .pred = Model_piCurve(fit$par, fit$model, data = df)
#' )
#' plot(df_ci$irradiance, df_ci$.obs, xlab = "Irradiance", ylab = "Photosynthesis Rate",
#'      ylim = c(min(df_ci$ymin), max(df_ci$ymax)))
#' lines(df_ci$irradiance, df_ci$.pred)
#' lines(df_ci$irradiance, df_ci$ymin, lty = "dashed")
#' lines(df_ci$irradiance, df_ci$ymax, lty = "dashed")
#'
addCI_to_piPred <- function(Fitted_Model, irrad, n_cores = 2, n_sample = 1000, Clevel = 0.95) {

    if (!("StDev" %in% names(Fitted_Model$par))){
        stop("No error estimates found. Fit model with `STATapp = 'MLE'` using Fit_piModel().")
    }

    # Get confidence intervals for each parameter
    CI_pars <- ConfInt_piCurve(Fitted_Model = Fitted_Model, Clevel = Clevel)

    if (length(which(is.na(CI_pars))) > 0){
        stop("NA detected in CI output - likely from numerical issues in solving the information matrix. See InfoMat_piCurve().")
    }
    # Filter parameter names
    par_name <- names(Fitted_Model$par)
    par_name <- par_name[par_name != "StDev"]

    # Extract CI bounds for each parameter into a named list
    CI_intervals <- setNames(
        lapply(par_name, function(p) CI_pars[paste0(c("L_CI.", "U_CI."), p)]),
        paste0(par_name, "_intv")
    )

    # Initialize parameter samples using uniform sampling from CI bounds
    set.seed(123)
    param_samples <- tibble(
        Pmax = runif(n_sample, CI_intervals$Pmax_intv[[1]], CI_intervals$Pmax_intv[[2]]),
        alpha = runif(n_sample, CI_intervals$alpha_intv[[1]], CI_intervals$alpha_intv[[2]])
    )

    # Sample optional parameters if available
    if ("beta" %in% par_name) {
        beta_intv <- CI_intervals$beta_intv
        param_samples$beta <- runif(n_sample, beta_intv[[1]], beta_intv[[2]])
    }

    if ("R" %in% par_name) {
        R_intv <- CI_intervals$R_intv
        param_samples$R <- runif(n_sample, R_intv[[1]], R_intv[[2]])
    }

    if ("shape" %in% par_name) {
        shape_intv <- CI_intervals$shape_intv
        param_samples$shape <- runif(n_sample, shape_intv[[1]], shape_intv[[2]])
    }

    # 2 tailed confidence interval quantile
    TTtest <- 1 - (1 - Clevel) / 2

    # Predict P for each sampled parameter set at each irradiance
    ci_summary <-
        mclapply(1:n_sample, function(i) {
            sampled_par <- setNames(as.numeric(param_samples[i,]), names(param_samples))

            pred_vals <- Model_piCurve(
                parameters = sampled_par,
                model_name = Fitted_Model$model,
                data = tibble(I = irrad)
            )

            return(tibble(irradiance = irrad, pred = pred_vals, draw = i))
        },
        mc.cores = n_cores) |>
        bind_rows() |>
        group_by(irradiance) |>
        summarise(
            ymin = quantile(pred, (1 - TTtest)),
            ymax = quantile(pred, TTtest),
            .groups = "drop"
        )

    return(ci_summary)
}
