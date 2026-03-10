#' Gregg and Carder (1990) atmospheric coefficients
#'
#' Spectral atmospheric coefficients and extraterrestrial solar irradiance
#' values used in the Gregg and Carder (1990) clear-sky maritime irradiance
#' model.
#'
#' @format A data frame with 351 rows and 5 variables:
#' \describe{
#'   \item{lambda}{Wavelength in vacuo (nm).}
#'   \item{H_o}{Spectral mean extraterrestrial solar irradiance
#'   (W m\eqn{^{-2}} nm\eqn{^{-1}}).}
#'   \item{a_oz}{Ozone absorption coefficient (cm\eqn{^{-1}}).}
#'   \item{a_w}{Water vapor absorption coefficient (cm\eqn{^{-1}}).}
#'   \item{a_o}{Oxygen absorption coefficient (cm\eqn{^{-1}}).}
#' }
#'
#' @details
#' This dataset contains wavelength-specific coefficients used by the
#' Gregg and Carder (1990) parameterization for estimating clear-sky
#' spectral irradiance at the sea surface.
#'
#' @source
#' Gregg, W. W., & Carder, K. L. (1990). A simple spectral solar irradiance
#' model for cloudless maritime atmospheres. \emph{Journal of Geophysical Research}.
#'
#' Source values were transcribed into CSV format for use in `piCurve`.
#'
#' @examples
#' head(GreggCarder1990)
"GreggCarder1990"
