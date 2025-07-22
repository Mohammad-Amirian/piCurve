#' Suggest Initial Parameter Values for P–I Curve Fitting
#'
#' Provides two sets of initial parameter estimates (`Pmax`, `alpha`, and `beta`) for fitting
#' photosynthesis–irradiance (PI) models. These values are computed based on the raw data
#' and are useful for initializing nonlinear curve fitting procedures.
#'
#' @param df A data frame containing columns `P` (photosynthesis rate) and `I` (irradiance).
#'
#' @return A matrix with two rows (`par1`, `par2`) and three columns (`Pmax`, `alpha`, `beta`)
#' representing candidate starting values for model parameters.
#'
#' @details
#' The function identifies the maximum photosynthesis rate (`Pmax`) and calculates two estimates
#' of the photochemical efficiency (`alpha`) and the photoinhibition term (`beta`) based on light
#' levels where photosynthesis exceeds 90% of `Pmax`. This heuristic approach helps provide
#' stable starting points for model optimization routines.
#'
#' @examples
#' params <- c(Pmax = 20, alpha = 0.6, beta = 0.1, R = 0)
#' df <- tibble::tibble(I = seq(0, 300, length = 50))
#' df$P <- Model_piCurve(parameters = c(params, beta = 0.3),
#'                       model_name = "2tanh",
#'                       data = df)
#' get_start_piPars(df)
#'
#'
#' @export
get_start_piPars <- function(df) {
    if (!all(c("P", "I") %in% names(df))) {
        stop("Input data frame must contain 'P' and 'I' columns.")
    }

    int_Pmax <- max(df$P, na.rm = TRUE)
    df_cutoff <- dplyr::filter(df, P / int_Pmax > 0.9)

    Imin <- min(df_cutoff$I, na.rm = TRUE)
    Imax <- max(df_cutoff$I, na.rm = TRUE)

    init_alpha <- int_Pmax / Imin
    init_beta <- int_Pmax / Imax

    output <- rbind(
        par1 = c(Pmax = int_Pmax, alpha = init_alpha, beta = init_beta),
        par2 = c(Pmax = int_Pmax, alpha = init_alpha, beta = 0)
    )

    return(output)
}
