
#' Confidence intervals
#'
#' @description
#' Compute standard errors and confidence intervals for the estimated parameter values for a provided PI model.
#'
#' @param Fitted_Model List -- the fitted model (output of \code{\link{Fit_piModel}} function).
#' @param Clevel Numeric -- the confidence level. Default 95%.
#'
#' @return Vector -- containing estimated errors along with a (\verb{Clevel})%
#' confidence interval for the estimated parameters. Default is 95%.
#' @export
#'
#' @importFrom stats qnorm
#' @examples
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = c(Pmax = 20, alpha = 0.6),
#'                  model_name = "LS5", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(STATapp = "MLE", Hessian = TRUE, data = df)
#'
#' # Calculate 95 % CI for the estimated parameters
#' ConfInt_piCurve(Fitted_Model = fit, Clevel = 0.95)
#'
ConfInt_piCurve <- function(Fitted_Model, Clevel = 0.95){
    n <- length(Fitted_Model$par)

    Est.Error <-  # solve the calculated Hessian matrix
        tryCatch(sqrt(diag(solve(Fitted_Model$hessian))),
                 error = function (e) rep(NA, n))

    # in case if the solution of the above matrix becomes NA
    names(Est.Error) <- names(Fitted_Model$par)

    # Two-tailed test; e.g for Clevel = 95%, TTtest = 1 - (1-0.95)/2
    TTtest <- 1 - (1-Clevel)/2
    Zscore <- stats::qnorm(TTtest)

    L_CI_val <- Fitted_Model$par - Zscore*Est.Error       # lower CI
    U_CI_val <- Fitted_Model$par + Zscore*Est.Error       # upper CI

    return(c(Est_Err = Est.Error, L_CI = L_CI_val, U_CI = U_CI_val))
}



#' Confidence intervals
#'
#' @description
#' Compute standard errors and confidence intervals for the estimated parameter values
#' when \verb{Hessian = TRUE} is forgotten in \code{\link{Fit_piModel}} function.
#'
#' @param Estimated_Params Vector -- the (optimal) parameter of a given PI model.
#' Further details \code{\link{Fit_piModel}}
#' @param model_name String – which model? (List of available models is given in \code{\link{Model_piCurve}})
#' @param data Vector – containing the irradiance profile.
#' @param Clevel Numeric -- the confidence level. Default 95%.
#'
#' @return Vector -- containing estimated errors along with a (\verb{Clevel})%
#' confidence interval for the estimated parameters. Default is 95%.
#' @export
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0, StDev = 2)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "LS5", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(params, model_name = "LS5", STATapp = "MLE", data = df)
#'
#' # Calculate CI for the fit
#' ReCal_CI_piCurve(Estimated_Params = fit$par, model_name = fit$model, data = df)
#'
#' @importFrom stats qnorm
ReCal_CI_piCurve <- function(Estimated_Params,
                             model_name,
                             data,
                             # data_type = c("light-limited", "light-saturating", "photoinhibition"),
                             Clevel = 0.95) {
    InfoMatrix <- InfoMat_piCurve(Estimated_Params, model_name, data)

    model_info  = c(par = list(Estimated_Params),
                    hessian = list(InfoMatrix))

    ConfInt_piCurve(Fitted_Model = model_info, Clevel)

}
