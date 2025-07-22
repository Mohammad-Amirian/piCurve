#####################################################################
# Statistical Functions:: MSE, adjusted_R2, AICc, BIC
#####################################################################

# MSE function ----
#' Statistical Adequacy Test
#'
#' @description
#' Wrapper to calculate Mean Squared Error (MSE) value.
#'
#' @param data Data frame –- containing the true photosynthesis rate profile.
#' @param model_fit Vector – containing the predicted photosynthesis rate profile.
#'
#' @return \samp{MSE} calculates measures the average squared difference between
#' predicted and actual values.
#'
#' @export
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(data = df)
#'
#' # Calculate the predicted photosynthetic rate
#' Phat <- Model_piCurve(parameters = fit$par, model_name = fit$model, data = df)
#'
#' MSE_piCurve(df, Phat)
#'
MSE_piCurve <- function(data, model_fit){
    sum( (data$P - model_fit)^2 / length(model_fit) )
}


# R2 and adjusted R2 ----
#' Statistical Adequacy Test
#'
#' @description
#' Wrapper to calculate R2 and adjusted R2 value.
#'
#' @param data Data frame –- containing the true photosynthesis rate profile.
#' @param model_fit Vector – containing the predicted photosynthesis rate profile.
#' @param Nparams Numeric -- total number of parameters used in the model.
#'
#' @return The function returns a vector consisting of R2 and adjusted R2.
#' @export
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(data = df)
#'
#' # Calculate the photosynthetic rate profile using the optimal values
#' model_fit <- Model_piCurve(parameters = fit$par, model_name = fit$model, data = df)
#'
#' R2_piCurve(df, model_fit, Nparams = length(fit$par))
#'
R2_piCurve <- function(data, model_fit, Nparams){

    # empirical photosynthetic rate profile
    P <- data$P
    # square difference between empirical P and estimated P
    diff_squared_val <- (P - model_fit)^2
    # take column-wise mean from the empirical data
    mean_databar <- mean(P, na.rm = TRUE)
    # subtract each vector from its mean value
    SST <- (P - mean_databar)^2
    # take column-wise sum
    SST <- sum(SST)
    # take column-wise sum from square diff between empirical data and the fit
    SSE <-  sum(diff_squared_val, na.rm = TRUE)

    # calculate R2
    R2 <- 1 - SSE/SST

    # calculate adjusted R2
    N = length(P)                # number of data points

    R2adj <- 1 - ( (N-1) / (N-Nparams) ) * (1-R2)

    return(c(R2 = R2, R2adj = R2adj))
}

# Calculate AIC, AICc, and BIC ----
#' Statistical Adequacy Test
#'
#' @description
#' Wrapper to calculate AIC, AICc, and BIC values.
#'
#' @param Fitted_Model List -- the fitted model (output of \code{\link{Fit_piModel}} function).
#' @param SampleN Numeric -- length of sample size.
#'
#' @return The function returns a vector that consists of AIC, AICc and BIC values.
#' @export
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(parameters = c(params, StDev = 2), model_name = "tanh",
#'                       STATapp = "MLE", Hessian = TRUE, data = df)
#'
#' # Calculate 95 % CI for the estimated parameters
#' AIC_AICc_BIC_piCurve(Fitted_Model = fit, SampleN = length(df$P))

AIC_AICc_BIC_piCurve <- function(Fitted_Model, SampleN){

    p = length(Fitted_Model$par)

    # AIC = 2(p - ln(likelihood)). To maximize the likelihood value, logL is
    # mutiplied by (-1) in MLE_fn, embedded in OptPIparams function. The negative
    # sign has thus been considered.
    AIC  <- 2 * (p + Fitted_Model$value)

    AICc <- AIC + 2 * p * (p + 1) / (SampleN - p -1)

    # BIC = p ln(SampleN) - 2 ln(likelihood)). Same as AIC applies to BIC
    BIC <- p * log(SampleN) + 2 * Fitted_Model$value

    return(c(AIC = AIC, AICc = AICc, BIC = BIC))
}


