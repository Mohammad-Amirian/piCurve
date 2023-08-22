
#' Model parameter confidence intervals
#'
#' @param Fitted_Model List -- the fitted model (output of \verb{OptPIparams} function).
#' @param Clevel Numeric -- the confidence level. Default 95%.
#'
#' @return Vector -- containing estimated errors along with a (\verb{Clevel})% confidence interval for the estimated parameters. Default is 95%.
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
#' df$PP <- # generate PP using Baly's rectangular hyperbola model
#'    Model_piCurve(parameters = params, model_name = "Eq3-Baly", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add some noise to PP
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- OptPIparams(parameters = c(params, StDev = 2), model_name = "Eq3-Baly",
#'                    STATapp = "MLE", Hessian = TRUE, data = df)
#'
#' # Calculate 95 % CI for the estimated parameters
#' ConfInt_piModel(Fitted_Model = fit, Clevel = 0.95)
#'
ConfInt_piModel <- function(Fitted_Model, Clevel = 0.95){
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


