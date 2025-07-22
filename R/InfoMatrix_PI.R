
#' Estimate the information matrix for PI models
#' @description
#' A wrapper to calculate the information matrix for PI models based on \code{\link{optimHess}} function.
#'
#' @param parameters Vector -- containing the optimal values listed below:
#' \itemize{
#'      \item{\code{Pmax} \eqn{\hspace{0.1cm}}: }{Maximum photosynthesis rate normalised by Chl_a (sign +),}
#'      \item{\code{alpha}: }{Light-saturation slope at low light level (sign +),}
#'      \item{\code{beta}  \eqn{\hspace{0.15cm}}: }{Photoinhibition rate at high light level (sign +),}
#'      \item{\code{R}  \eqn{\hspace{0.7cm}}: }{dark reaction parameter (sign + OR -),}
#'      \item{\code{shape}  \eqn{\hspace{0cm}}: }{shape parameter (sign +) -- Only required for some of the models. See \samp{Details}}.
#'      \item{\code{StDev}  \eqn{\hspace{0cm}}: }{standard deviation parameter -- Only required when \samp{STATapp = MLE}}.
#' }
#' @param model_name String -- which model? (See \samp{details})
#'
#' @param data Vector -- Containing the irradiance profile.
#' @param GradientFn String -- A function that produces the gradient for \verb{BFGS},
#' @param ... Further arguments to be passed to fn and gr.
#' @param Control List -- a list of control parameters. Further details \code{\link{optim}}.
#'
#' @return A symmetric matrix that estimates the information/hessian matrix at the optimal solution.
#'
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
#'    Model_piCurve(parameters = params, model_name = "tanh", data = df) +
#'    5 * rnorm(25, 0, 0.25)    # add some noise to P
#'
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' fit <- Fit_piModel(params, model_name = "tanh", STATapp = "MLE", data = df)
#'
#' # Calculating the information matrix
#' InfoMat_piCurve(parameters = fit$par, model_name = fit$model, data = df)
#'
#' # different accuracy
#' InfoMat_piCurve(parameters = fit$par, model_name = fit$model, data = df,
#'                 Control = list( ndeps = rep(1e-4, length = length(fit$par)) ) )
#'
#'InfoMat_piCurve(parameters = fit$par, model_name = fit$model, data = df,
#'                Control = list( ndeps = c(rep(1e-4, 3), StDev = 1e-3)) )

#' @importFrom stats optimHess

InfoMat_piCurve <-
    function(parameters, model_name, data, GradientFn = NULL, ..., Control = list()) {
    # data format check ----
    if (!all(c("P", "I") %in% names(data))) {
        print("Input data frame must contain 'P' and 'I' columns. Let's fix the issue ...")
        data <- FormatCheck_piCurve(data = data)
    }

    # data type = ? ----
    data_type <-
        ifelse(
            nameFormat(model_name) %in% Pool_eqName_lm, "ll", # light-limited
            ifelse(nameFormat(model_name) %in% Pool_eqName_ls,
                   "ls", # light-saturating
                   "ph" # photoinhibition
            )
        )

    # select the model specified by the user
    equation <- Model_setup(which_model = model_name, data_type)


    Cal_logL <- function(parameters, data) {
        # predicted P for a given parameters
        Phat <- equation(parameters = parameters, data)

        # calculate log-likelihood function and return it
        logl <-
            LogL(data, model_fit = Phat, StDev = parameters["StDev"])

        # return negative logl as it gets maximized using optim
        return(-logl)
    }

    info_matrix <-
        optimHess(
            par = parameters,
            fn = Cal_logL,
            gr  = GradientFn,
            ...,
            control = Control,
            data = data
        )

    return(info_matrix)

}

