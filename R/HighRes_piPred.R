#' Generate High-Resolution Model Predictions for PI Curves
#'
#' @description
#' This function takes a fitted photosynthesis–irradiance (PI) model and returns
#' a high-resolution prediction of photosynthetic rate (`P`) across a finely spaced
#' irradiance (`I`) range based on the input dataset.
#'
#' @param Fitted_Model A fitted model object returned by \code{Fit_piModel()}, containing the parameter estimates and model name.
#' @param data A data frame or tibble with at least columns `I` (irradiance) and `P` (photosynthesis). Used to define the prediction range.
#' @param length_out Integer. Number of points to generate across the irradiance range. Default is 50.
#'
#' @return A tibble with high-resolution irradiance values and the corresponding model-predicted photosynthesis rates (`Phat`).
#'
#' @details
#' If the input `data` does not contain the expected columns `I` and `P`, the function attempts to format it using `FormatCheck_piCurve()`.
#'
#' @examples
#' df <- tibble::tibble(
#'   I = c(0, 33, 67, 100, 133, 167, 200),
#'   P = c(0, 15, 19, 20, 20, 20, 20)
#' )
#'
#' fit <- Fit_piModel(data = df)
#'
#' pred_highres <- highRes_piPred(fit, data = df)
#'
#' @export
highRes_piPred <- function(
        Fitted_Model,
        data,
        length_out = 50  # desired number of points in high-res output
) {
    # Check if input data has required columns ----
    if (!all(c("P", "I") %in% names(data))) {
        message("Input data frame must contain 'P' and 'I' columns. Attempting to format...")
        data <- FormatCheck_piCurve(data = data)
    }

    # Generate evenly spaced irradiance values within observed range ----
    df_high_res <- tibble(
        I = seq(
            from = min(data$I),
            to = max(data$I),
            length.out = length_out
        ),

        # Predict photosynthesis using the fitted model at each irradiance value
        Phat = Model_piCurve(
            parameters = Fitted_Model$par,
            model_name = Fitted_Model$model,
            data = tibble(I = I)  # provide irradiance values to prediction function
        )
    )

    return(df_high_res)
}
