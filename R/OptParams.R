
#' Fit PI models
#' @description
#' \samp{Fit_piModel} uses general-purpose optimization to identify the optimal parameters
#' in the photosynthesis-irradiance curve models.
#' The method includes the \samp{Nelder-Mead}, \samp{quasi-Newton}, and \samp{conjugate-gradient} algorithms,
#' as well as options for box-constrained optimization and simulated annealing. Two different
#' statistical approaches, Mean Squared Error (MSE) and Maximum Likelihood Estimation (MLE),
#' are available for objective function optimization.
#'
#' @param parameters Vector -- containing the values listed below. If `NULL`, default values are generated using \code{\link{get_start_piPars}}.
#' \itemize{
#'      \item{\code{Pmax} \eqn{\hspace{0.1cm}}: }{Maximum photosynthetic rate normalised by Chl_a (sign +),}
#'      \item{\code{alpha}: }{Light-saturation slope at low light level (sign +),}
#'      \item{\code{beta}  \eqn{\hspace{0.15cm}}: }{Photoinhibition parameter (sign +),}
#'      \item{\code{R}  \eqn{\hspace{0.7cm}}: }{Dark reaction rate (\samp{Optional} with sign + OR -)}
#'      \item{\code{shape}  \eqn{\hspace{0cm}}: }{Shape parameter (sign +) -- Only required for some of the models. See \samp{Details}}.
#'      \item{\code{StDev}  \eqn{\hspace{0cm}}: }{Standard deviation parameter -- Only required when \samp{STATapp = MLE}}.
#' }
#' @param model_name String -- which model? If `NULL`, defaults to `Ph10` (double-tanh) model with and without photoinhibition. See \code{\link{Model_piCurve}} and \samp{details}.
#' @param STATapp String -- Providing the statistical methods for optimization,
#' which are \samp{MSE} or \samp{MLE}. See \samp{Details}. Default is \samp{MSE}.
#' @param data Data frame -- Containing  `I` (irradiance) and `P` (photosynthesis rate) columns
#' @param GradientFn String -- A function that produces the gradient for \verb{BFGS},
#' \verb{CG} and \verb{L-BFGS-B} techniques, and if it is not provided, a numerical
#' estimation using \samp{finite-differences} will be applied.
#' @param ... An option to pass further arguments to \verb{GradientFn}.
#' @param Method String -- The method employed to determine the optimal values for the
#' model's parameters. Further details \code{\link{optim}}.
#' @param LowerBnd Numeric -- Lower boundaries on the variables for the \verb{L-BFGS-B}
#' method, or bounds in which to search for method \verb{Brent}.
#' @param UpperBnd Numeric -- Upper boundaries on the variables for the \verb{L-BFGS-B}
#' method, or bounds in which to search for method \verb{Brent}.
#' @param Control List -- a list of control parameters. Further details \code{\link{optim}}.
#' @param Hessian Logical -- Should the Hessian matrix be returned through numerical
#' differentiation? Default is \verb{FALSE}.
#'
#' @return Below is a list of the results obtained from the optimization.
#' \itemize{
#' \item{\verb{par:}} The optimal parameters obtained through optimization by
#' using either the mean squared error (\samp{MSE}) or maximum likelihood estimate (\samp{MLE}) method.
#'
#' \item{\verb{value:}}
#' The MSE value on the primary production (\verb{P}) profile when
#' \verb{STATapp = MSE} is met or \verb{STATapp} is not specified. This metric
#' measures the average squared difference between predicted and actual values
#' in regression analysis. A lower value for the mean MSE indicates a better fit.
#' For \verb{STATapp = MLE}, the optimization uses Maximum Likelihood Estimation
#' method to find the optimal parameters. The value in this case stands for the
#' likelihood of observing the primary production (\verb{P}) for the estimated
#' parameters. A higher value for the MLE indicates a better fit as the method
#'  aims to find the parameter values that make the observed data most
#' likely to occur according to the assumed distribution
#'
#' \item{\verb{counts:}}
#' This is a vector consisting of two integers that indicate the number of times
#' the model and its gradient \verb{GradientFn} were called during a certain process.
#' This count does not include the calls made to compute the Hessian (if requested),
#' or any calls to the model to compute a finite-difference approximation to the gradient.
#'
#' \item{\verb{message:}}
#' A character string that provides any extra information given by the optimizer
#' upon completion of a process. if there is no such additional information, a value of \verb{NULL} will be returned.
#'
#' \item{\verb{hessian:}}
#' A symmetric matrix that estimates the Hessian at the optimal solution, but only if the hessian argument is set to \verb{TRUE}.
#' It is important to note that this is the Hessian of the unconstrained problem, even if there are active box constraints.
#'
#' \item{\verb{convergence:}}
#' An integer value representing the optimization process is completed. A value of 0 indicates
#' that the process has completed successfully, which is always the case for \verb{SANN} and \verb{Brent}.
#' However, there may be other possible error codes listed below, which could be generated during the process.
#'
#' \itemize{
#'      \item{\verb{1 :}  }{indicates that the iteration limit maxit had been reached.}
#'      \item{\verb{10:}}{indicates that the \verb{Nelder-Mead} simplex method has encountered degeneracy, meaning that
#'      two or more vertices of the simplex have converged to the same point.}
#'      \item{\verb{51:}}{indicates a warning that has been issued by the  \verb{L-BFGS-B} method.
#'      Further details regarding this warning can be found in the \verb{message} component of the output.}
#'      \item{\verb{52:}} {indicates an error that has occurred within the \verb{L-BFGS-B} method.
#'      Additional details about this error can be found in the \verb{message} component of the output.}
#' }
#'
#' \item{\verb{info_criteria:}} A vector consisting AIC, AICc, and BIC values when \verb{MLE} approach is chosen.
#' The main difference between these values lie in the penalties they impose on model complexity.
#' AIC tends to select more complex models compared to BIC, while AICc provides a correction
#' to AIC for small sample sizes. BIC, on the other hand, penalizes model complexity more
#' heavily than AIC and AICc, often resulting in the selection of simpler models.
#'
#' \item{\verb{model:}} Name of equation used to model the given data sample.
#'
#' \item{\verb{SQA:}} Statistical Quality Assessment consists of a vector incorporating
#' R2 and adjusted R2 values. R2 values provide information about the
#' goodness of the fit by indicating the proportion of data variation that is
#' explained by the model. Adjusted R2 considers a penalty for over parameterization.
#' Together, these two measures help evaluating the accuracy
#' of the model's predictions and provide insights into how well the model
#' captures/explains the observed data variation. Note that these values should
#' not be used as a model selection criteria. AICc or BIC should be used for that
#' purpose.
#' }
#'
#' @details
#'
#' The \strong{\verb{STATapp}} argument refers to the statistical approach.  To achieve
#' a more accurate estimation in a regression model (greater R\eqn{^2}), users may consider
#' using the MSE method over the MLE method, as the former measures the average
#' squared difference between predicted and actual values. In contrast, the MLE
#' method determines the best possible values by maximizing the likelihood of
#' observing the data within the given model, thus offering the user a detailed
#' insight into the estimation error of each parameter. To get such insight when
#' \verb{STATapp = MLE}, the \verb{Hessian} argument should be set to \verb{TRUE}. If forgotten,
#' the user can use \verb{optimHess} function.
#'
#' This package includes 36 models which can be found in \code{\link{Model_piCurve}} both
#' commonly-used and recently-developed PI models. The list of models is provided below.
#'
#' This function utilizes the \verb{optim()} function from the \verb{stats}
#' package to perform the optimization. As such, the function makes available all
#' of the arguments that can be used with the \verb{optim()} function within itself.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 200, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'   Model_piCurve(parameters = params, model_name = "LS5", data = df) +
#'   5 * rnorm(25, 0, 0.25)    # add noise
#'
#' # Estimate the optimal values for the generated dataset using MSE method
#' # Example 1: manual initial values
#' Fit_piModel(parameters = c(Pmax = 19, alpha = 0.5, R = 0), model_name = "LS5", data = df)
#'
#' # Example 2: automated initial values and model selection
#' Fit_piModel(data = df)
#'
#' # Estimate the optimal values for the generated dataset using MLE method
#' Fit_piModel(STATapp = "MLE", Hessian = TRUE, data = df)
#'
#' @importFrom fBasics Heaviside
#' @importFrom stats optim
#'
Fit_piModel <- function(data,
                           parameters = NULL,
                           model_name = NULL,
                           STATapp = c("MSE", "MLE"),
                           GradientFn = NULL,
                           ...,
                           Method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                           LowerBnd = -Inf,
                           UpperBnd =  Inf,
                           Control = list(),
                           Hessian = FALSE) {

    # data format check ----
    if (!all(c("P", "I") %in% names(data))) {
        print("Input data frame must contain 'P' and 'I' columns. Let's fix the issue ...")
        data <- FormatCheck_piCurve(data = data)
    }

    # dataType ? ----
    if(is.null(model_name)){
        data_type <- DataType_piCurve(data = data, n_cores = 1)$data_type

        model_name <-
            ifelse(data_type == "ph", "ph10",
                   ifelse(data_type == "ls", "ls5", "ll"))
    } else {
        model_name <- tolower(model_name)
        data_type <-
            ifelse(
                nameFormat(model_name) %in% Pool_eqName_lm,
                "ll", # light-limited
                ifelse(
                    nameFormat(model_name) %in% Pool_eqName_ls,
                    "ls", # light-saturating
                    "ph" # photoinhibition
                )
            )
    }

    # over-write!
    if (data_type == "lm") {
        fit_lm <- summary(lm(I ~ P, data))
        stop(paste("Light-limited data detected. Model summary:\n",
                   paste(capture.output(fit_lm), collapse = "\n")))
    }

    # if parameters is missing, use get_start_piPars() function ----
    if (is.null(parameters)) {
        parameters <- cbind(get_start_piPars(data), StDev = 1, shape = 0.9)[1,]
    }

    # If STATapp is "MLE" and StDev is missing from parameters, set default StDev = 1
    ifelse(
        missing(STATapp) || tolower(STATapp) == "mse",
        parameters,
        ifelse(
            "StDev" %in% names(parameters),
            parameters,
            parameters <- c(parameters, StDev = 1)
        )
    )

    # select the model specified by the user
    equation <- Model_setup(which_model = model_name, data_type)

    # MSE_fn needs an equation name, requiring us to define the function in
    # the global environment due to Scoping Rules of R
    MSE_fn <- function(parameters, data) {
        Phat <- # predicted P for a given parameters
            tryCatch(
                equation(parameters = parameters, data),
                error = function (e)
                    rep(0, data$P)
                )

        MSE_piCurve(data, model_fit = Phat) # return MSE value
    }


    # same comment as MSE_n concerning scoping rules of R applies to MLE_fn
    MLE_fn <- function(parameters, data) {
        # predicted P for a given parameters
        Phat <-
            tryCatch(
                equation(parameters = parameters, data),
                error = function (e)
                    rep(0, data$P)
                )

        logl <- # calculate log-likelihood function
            tryCatch(
                LogL(data, model_fit = Phat, StDev = parameters["StDev"]),
                error = function (e) -1e5  # generate a large number in case of err
                )

        # return negative logl as it gets maximized using optim
        return(-logl)
    }

    # find the optimal parameters using optim fn with either MSE or MLE approach.
    # Default choice MSE
    ifelse(
        missing(STATapp) || tolower(STATapp) == "mse",
        fit <- optim(
            par = parameters,
            fn = MSE_fn,
            data = data,
            gr  = GradientFn,
            ...,
            method  = Method,
            lower   = LowerBnd,
            upper   = UpperBnd,
            control = Control,
            hessian = Hessian
            ),
        ifelse(
            tolower(STATapp) == "mle",
            fit <- optim(
                par = parameters,
                fn = MLE_fn,
                data = data,
                gr  = GradientFn,
                ...,
                method  = Method,
                lower   = LowerBnd,
                upper   = UpperBnd,
                control = Control,
                hessian = Hessian
                )
            ,
            stop(
                "Please ensure that the STATapp field is either empty or set to 'mse' or 'mle'"
            )
        )
    )

    # in case if any letter exists in model name (provided by user), remove them
    # model_name <- gsub("[^1,2,3,4,5,6,7,8,9,10]", "", model_name)

    # Other than dark resperation rate (R), rest of the params are pushed into abs() in the main
    # functions, making negative and positive values equivalent. There is two
    # exceptions associates with LS6_Prioul and Ph08_Prioul in which the shape parameter belongs
    # to (-Inf, Inf). Thus, I am excluding them from other models before taking
    # the absolute value (abs()) of the optimized values.

    ifelse(nameFormat(model_name) %in% c("ls6", "ph08"),
           # In equations ls6 and ph08, shape ranges from -Inf to Inf
           fit$par[!(names(fit$par) %in% c("shape", "R"))] <- abs(fit$par[!(names(fit$par) %in% c("shape", "R"))]),
           # the other the equations
           fit$par[names(fit$par) != "R"] <- abs(fit$par[names(fit$par) != "R"])
           )

    ifelse(
        missing(STATapp) || tolower(STATapp) == "mse",
        fit,
        # add AIC, AICc, and BIC info to the fit for MLE method
        fit$info_criteria <-
            AIC_AICc_BIC_piCurve(Fitted_Model = fit, SampleN = length(data$P))
    )

    # --- clean up par depending on the method ----
    ifelse(
        missing(STATapp) || tolower(STATapp) == "mse",
        fit$par <- fit$par[names(fit$par) != "StDev"], # rm StDev
        fit
    )

    ifelse(
        data_type == "ls",
        fit$par <- fit$par[names(fit$par) != "beta"], # rm beta
        ifelse(
            data_type == "ll",
            fit$par <- fit$par[!(names(fit$par) %in% c("beta", "alpha"))],
            fit
        )
    )

    ### ----
    if (model_name %in% steel_model_names) {
        fit$par <- fit$par[names(fit$par) != "beta"]
    }

    if (!(model_name %in% extra_param_models)) {
        fit$par <- fit$par[names(fit$par) != "shape"] # rm shape param
    }

    # update hessian matrix
    if ("StDev" %in% names(fit$par)){
        idx = which(colnames(fit$hessian) %in% names(fit$par))
        fit$hessian <- fit$hessian[idx, idx]
    }
    # ----

    # add the model name to the optimization results
    fit$model <- model_name

    # add R2 ad adjusted R2 values to the optimization results
    Phat <- equation(parameters = fit$par, data)
    fit$SQA <- R2_piCurve(data, model_fit = Phat, Nparams = length(fit$par))

    return(fit)
}


