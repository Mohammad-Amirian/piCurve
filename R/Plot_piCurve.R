#' Plot Fitted PI Curve with Optional Confidence Intervals
#'
#' Visualizes a fitted photosynthesis–irradiance (PI) model along with optional confidence intervals on
#' observed data points.
#'
#' @param Fitted_Model A model object returned by \code{\link{Fit_piModel}} containing parameter estimates and model type.
#' @param data A data frame containing observed data with columns \code{I} (irradiance) and \code{P} (photosynthesis).
#' @param add_CI Logical. Whether to add confidence intervals around the fitted curve (default: \code{FALSE}, see details).
#' @param n_cores Integer. Number of CPU cores to use if computing confidence intervals (default: 1).
#' @param length_out Integer. Number of points to use when generating the high-resolution fitted curve (default: 50).
#' @param point_size Numeric. Size of observed data points in the plot (default: 2).
#'
#' @return A \code{ggplot} object showing the observed data points, fitted PI curve, and optionally the confidence interval bounds.
#'
#' @details
#' If `add_CI = TRUE`, the function internally calls \code{\link{addCI_to_piPred}} to compute confidence intervals at each irradiance level used for prediction.
#'
#' If the input data is not properly formatted (missing columns "P" or "I"), the function attempts to fix it using \code{\link{FormatCheck_piCurve}}.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_bw
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   I = c(0, 33, 67, 100, 133, 167, 200),
#'   P = c(0, 15, 19, 20, 20, 20, 20)
#' )
#'
#' fit <- Fit_piModel(
#'   parameters = c(Pmax = 19, alpha = 0.5, R = 0, StDev = 0.1),
#'   data = df, STATapp = "MLE", Hessian = TRUE
#' )
#'
#' Plot_piCurve(fit, data = df, add_CI = TRUE, n_cores = 2)
#'
Plot_piCurve <- function(
        Fitted_Model,
        data,
        add_CI = FALSE,
        n_cores = 1,
        length_out = 50,
        point_size = 2
) {
    # Ensure data has proper structure
    if (!all(c("P", "I") %in% names(data))) {
        message("Input data frame must contain 'P' and 'I' columns. Attempting to fix...")
        data <- FormatCheck_piCurve(data = data)
    }

    # High-resolution prediction
    df_high_res <- highRes_piPred(Fitted_Model, data, length_out)


    # Add confidence interval if requested
    if (add_CI == TRUE) {
        # suppress warnings inside this call
        ci <- suppressWarnings(ConfInt_piCurve(Fitted_Model))
        if (any(is.na(ci))) {
            df_high_res <- df_high_res

            message("CI ignored: NA detected while inverting the information matrix; see ConfInt_piCurve().")

        } else{

            df_ci <- addCI_to_piPred(Fitted_Model, irrad = df_high_res$I, n_cores = n_cores)
            df_high_res <- cbind(df_high_res, df_ci[, c("ymin", "ymax")])

        }
    }

    # Construct base plot
    pp <- ggplot2::ggplot(data, aes(x = I, y = P)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::geom_line(data = df_high_res, mapping = aes(x = I, y = Phat), col = "black", alpha = 0.9) +
        ggplot2::labs(x = "Irradiance", y = "Photosynthesis Rate") +
        ggplot2::theme_bw()

    # Add CI bounds if available
    if (add_CI == TRUE) {
        if (!any(is.na(ci))) {
            pp <- pp +
                ggplot2::geom_line(data = df_high_res, aes(x = I, y = ymin), linetype = 2, col = "black", alpha = 0.6) +
                ggplot2::geom_line(data = df_high_res, aes(x = I, y = ymax), linetype = 2, col = "black", alpha = 0.6)
        }
    }

    return(pp)
}
