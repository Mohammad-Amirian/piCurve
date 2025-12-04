#' Compute phytoplankton absorption coefficient
#'
#' Estimate chlorophyll-a-specific phytoplankton absorption \eqn{a^*_{ph}(\lambda)}
#' and phytoplankton absorption \eqn{a_{ph}(\lambda)} from chlorophyll-a concentration
#' using the Bricaud et al. (1995) parameterization.
#'
#' @param chla Numeric vector. Chlorophyll-a concentration in mg m\eqn{ ^{-3} }. Must be > 0.
#' @param coef_df Optional data frame providing columns `lambda_nm`, `A`, `B`.
#'   Defaults to [piCurve::bricaud1995_coef].
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{chla}{Chlorophyll-a (mg m\eqn{^{-3}}).}
#'   \item{lambda_nm}{Wavelength (nm).}
#'   \item{aStar_ph}{Chl-specific phytoplankton absorption, \eqn{a^*_{ph}(\lambda)}; unit: m\eqn{^2} mg\eqn{^{-1}}.}
#'   \item{a_ph}{Phytoplankton absorption, \eqn{a_{ph}(\lambda)}; unit: m\eqn{^{-1}}.}
#' }
#' @details
#' Parameterization:
#' \deqn{a^*_{ph}(\lambda) = A(\lambda)\ \mathrm{Chl}^{-B(\lambda)}}
#' \deqn{a_{ph}(\lambda) = a^*_{ph}(\lambda)\ \mathrm{Chl}}
#'
#' @references
#' Bricaud, A., Babin, M., Morel, A., & Claustre, H. (1995). Variability in the chlorophyll-specific
#' absorption coefficients of natural phytoplankton: Analysis and parameterization.
#' \emph{Journal of Geophysical Research: Oceans}, 100(C7), 13321–13332.
#'
#' @examples
#' absCoef(chla = 0.1)
#' absCoef(chla = c(0.05, 0.1, 1))
#'
#' @export
absCoef <- function(chla, coef_df = piCurve::bricaud1995_coef) {
    if (!is.numeric(chla)) {
        stop("`chla` must be numeric (mg m^-3).", call. = FALSE)
    }
    if (any(!is.finite(chla))) {
        stop("`chla` must contain only finite values.", call. = FALSE)
    }
    if (any(chla <= 0)) {
        stop("`chla` must be strictly > 0 (mg m^-3).", call. = FALSE)
    }

    req <- c("lambda_nm", "A", "B")
    if (!all(req %in% names(coef_df))) {
        stop("`coef_df` must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
    }

    # Avoid importing dplyr if you prefer base; but tibble is nice for return shape.
    # Use tidyr::crossing-like behavior via base expand.grid for minimal deps:
    out <- merge(
        data.frame(chla = as.numeric(chla)),
        data.frame(
            lambda_nm = as.numeric(coef_df$lambda_nm),
            A = as.numeric(coef_df$A),
            B = as.numeric(coef_df$B)
        ),
        by = NULL
    )

    out$aStar_ph <- out$A * (out$chla)^(-out$B)
    out$a_ph <- out$aStar_ph * out$chla

    # Return tidy tibble
    tibble::as_tibble(out[, c("chla", "lambda_nm", "aStar_ph", "a_ph")])
}
