
#' Confidence intervals
#'
#' @description
#' A wrapper to compute confidence intervals for the optimal values, in addition
#' to calculating the estimation errors, for a provided PI model.
#'
#' @param Fitted_Model List -- the fitted model (output of \verb{OptPar_piCurve} function).
#' @param Clevel Numeric -- the confidence level. Default 95%.
#'
#' @return Vector -- containing estimated errors along with a (\verb{Clevel})%
#' confidence interval for the estimated parameters. Default is 95%.
#' @export
#'
#' @importFrom stats qnorm
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$PP <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df, data_type = "ls") +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- OptPar_piCurve(parameters = c(params, StDev = 2), model_name = "tanh",
#'                      STATapp = "MLE", Hessian = TRUE, data = df, data_type = "ls")
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
#' A wrapper to compute confidence intervals for the optimal values, in addition
#' to calculating the estimation errors when \verb{Hessian = TRUE} is forgotten
#' in \verb{OptPar_piCurve} function.
#'
#' @param Estimated_Params Vector -- the (optimal) parameter of a given PI model.
#' Further details \code{?piCurve::OptPar_piCurve}
#' @param model_name String – which model? (List of available models is given in \code{?piCurve::Model_piCurve})
#' @param data Vector – containing the irradiance profile.
#' @param data_type String – data sample type?
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
#' df$PP <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df, data_type = "ls") +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- OptPar_piCurve(params, model_name = "tanh", STATapp = "MLE", data = df, data_type = "ls")
#'
#' # Calculate CI for the fit
#' ReCal_CI_piCurve(Estimated_Params = fit$par, model_name = fit$model, data = df, data_type = "ls")
#'
#' @importFrom stats qnorm
ReCal_CI_piCurve <- function(Estimated_Params,
                             model_name,
                             data,
                             data_type = c("light-limited", "light-saturating", "photoinhibition"),
                             Clevel = 0.95) {
    InfoMatrix <- InfoMat_piCurve(Estimated_Params, model_name, data, data_type)

    model_info  = c(par = list(Estimated_Params),
                    hessian = list(InfoMatrix))

    ConfInt_piCurve(Fitted_Model = model_info, Clevel)

}
