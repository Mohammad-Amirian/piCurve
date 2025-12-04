#' Phytoplankton absorption coefficients
#'
#' Coefficients used to estimate the chlorophyll-a-specific phytoplankton absorption spectrum
#' \eqn{a^*_{ph}(\lambda)} following Bricaud et al. (1995). The parameterization is typically
#' applied as:
#' \deqn{a^*_{ph}(\lambda) = A(\lambda)\ \mathrm{Chl}^{-B(\lambda)}}
#' and spectral phytoplankton absorption is:
#' \deqn{a_{ph}(\lambda) = a^*_{ph}(\lambda)\ \mathrm{Chl}}
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{lambda_nm}{Wavelength (nm).}
#'   \item{A}{Wavelength-specific coefficient (units consistent with \eqn{a^*_{ph}} in m\eqn{^2} mg\eqn{^{-1} }).}
#'   \item{B}{Wavelength-specific exponent (dimensionless).}
#'   \item{r2}{Coefficient of determination (\eqn{R^2}) of the log--log regression
#'   used to derive \eqn{A(\lambda)} and \eqn{B(\lambda)} at each wavelength.}
#' }
#'
#' @details
#' The dataset is intended for spectral calculations. If you need band-averaged or PAR-weighted
#' absorption, compute \eqn{a_{ph}(\lambda)} and then apply your weighting function externally.
#'
#' @source
#' Bricaud, A., Babin, M., Morel, A., & Claustre, H. (1995). Variability in the chlorophyll-specific
#' absorption coefficients of natural phytoplankton: Analysis and parameterization.
#' \emph{Journal of Geophysical Research: Oceans}, 100(C7), 13321–13332.
#'
#' @keywords datasets
"bricaud1995_coef"
