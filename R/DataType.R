#' Automatic Data Classifier
#'
#' Fits multiple photosynthesis–irradiance (PI) models to the data in parallel and infers the data type based on the best adjusted R-squared.
#'
#' @param data A data frame containing `I` (irradiance) and `P` (photosynthesis rate) columns.
#' @param n_cores Integer. Number of CPU cores to use for parallel model fitting. Default is 2 to comply with CRAN check limits.
#' @return A list with a single named element `data_type`, indicating the data type: _ll (light-limited)_, _ls (light-saturated)_, and _ph (photoinhibition)_.
#'
#' @details
#' The function fits a set of candidate models: photoinhibition models ("Ph09", "Ph10"),
#' light-saturating models ("LS1", "LS5"), and a linear model for light-limited cases
#' (models listed in \code{\link{Model_piCurve}}) and labels the data based on adjusted R2.
#' R²adj is chosen over other statstical metric, as the PI data could be extremely imbalanced (see \samp{References} for details).
#'
#' If the linear model fits exceptionally well (R²adj > 0.94), it is strongly favored.
#' For details on these model choices, see \samp{References}.
#'
#' @importFrom dplyr bind_rows arrange slice
#' @importFrom parallel mclapply detectCores
#' @importFrom stats dnorm lm na.omit
#' @importFrom utils capture.output
#' @export
#'
#' @references
#' Amirian M.M., V Finkel Z., Devred E., Irwin A.J.,
#' “\emph{A new parameterization of photoinhibition for phytoplankton},
#' arXiv (2024) 1–33. 10.48550/arXiv.2412.17923.
#'
#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, beta = 0.1)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 300, length = 25))
#'
#' df$P <- # generate the photosynthetic rate using Jassby-tanh (LS5) model
#'  Model_piCurve(parameters = params, model_name = "2tanh", data = df) +
#'  2 * rnorm(25, 0, 0.25)
#'
#'  DataType_piCurve(data = df)
#'
DataType_piCurve <- function(data, n_cores = 2) {

    data <- FormatCheck_piCurve(data)
    params <- get_start_piPars(data)
    list_models <- c("lm", "ls1", "ls5", "ph09", "ph10")

    df_out <- vector("list", nrow(params))

    for (x in seq_len(nrow(params))) {
        df_out[[x]] <-
            parallel::mclapply(
                seq_along(list_models),
                function(idx_model){
                    model <- list_models[idx_model]

                    param_vec <- params[x, ]

                    fit <- tryCatch(
                        Fit_piModel(
                            parameters = param_vec,
                            model_name = model,
                            data = data
                        ),
                        error = function(e)
                            NA
                    )


            # Evaluate adjusted R-squared
            R2adj <- 0  # default fallback

            if (!is.null(fit) && !all(is.na(fit))) {
                if (model %in% c("ph09", "ph10")) {
                    Ib <- as.numeric(fit$par["Pmax"] / fit$par["beta"])
                    R2adj <- ifelse(Ib > max(data$I, na.rm = TRUE), 0, fit$SQA["R2adj"])
                } else {
                    R2adj <- fit$SQA["R2adj"]
                }

                # over-write!
                if (model == "lm"){
                    fit_lm <- summary(lm(I ~ P, data))
                    R2adj <- fit_lm$adj.r.squared
                }

                if (model == "lm" && R2adj > 0.94) {
                    R2adj <- 1  # strongly favor pure linear if extremely good fit
                }
            }

            tibble::tibble(model = model, R2adj = as.numeric(R2adj))
        },
        mc.cores = n_cores
        )
    }

    # Flatten and find the best-fitting model
    df_flat <- bind_rows(unlist(df_out, recursive = FALSE))
    df_best <- df_flat |> arrange(model) |> slice(which.max(R2adj))

    return(list(data_type = gsub("[0-9]+", "", df_best$model) ))
}


