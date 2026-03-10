#' Calculate atmospheric optical mass terms
#'
#' Computes relative optical path lengths for air, pressure-corrected air,
#' and ozone as functions of solar zenith angle.
#'
#' @param zenith_angle Numeric scalar. Solar zenith angle in degrees.
#'   Must satisfy `0 <= zenith_angle < 90`.
#' @param non_stnd_pressure Numeric scalar. Atmospheric pressure in millibars (mb).
#'   Default is `1010.5`.
#'
#' @return A named numeric vector with:
#' \describe{
#'   \item{M_theta}{Relative air mass.}
#'   \item{M_theta_prime}{Pressure-corrected air mass.}
#'   \item{M_ozone}{Relative ozone optical mass.}
#' }
#'
#' @details
#' This function computes atmospheric mass terms used in the Gregg&Carder 1990
#' clear-sky spectral irradiance framework. `M_theta_prime` scales air mass by
#' the ratio of local pressure to standard pressure (`1013.25` mb).
#'
#' @examples
#' atmospheric_mass_piCurve(zenith_angle = 30)
#'
#' @export
atmospheric_mass_piCurve <- function(zenith_angle, non_stnd_pressure = 1010.5) {
    if (!is.numeric(zenith_angle) || length(zenith_angle) != 1L || is.na(zenith_angle)) {
        stop("`zenith_angle` must be a single non-missing numeric value.")
    }
    if (zenith_angle < 0 || zenith_angle >= 90) {
        stop("`zenith_angle` must be in [0, 90).")
    }
    if (!is.numeric(non_stnd_pressure) || length(non_stnd_pressure) != 1L || is.na(non_stnd_pressure)) {
        stop("`non_stnd_pressure` must be a single non-missing numeric value.")
    }

    # air mass (approximation for the path length through the atmosphere)
    theta_rad <- zenith_angle * pi / 180
    cos_zenith <- cos(theta_rad)
    M_theta <- 1 / (cos_zenith + 0.15 * ((93.885 - zenith_angle)^(-1.253)))

    P0 <- 1013.25

    # corrected air mass
    M_theta_prime <- M_theta * non_stnd_pressure / P0

    # approximation for the path length through the ozone
    M_ozone <- 1.0035 / sqrt(cos_zenith^2 + 0.007)

    c(
        M_theta = M_theta,
        M_theta_prime = M_theta_prime,
        M_ozone = M_ozone
    )
}
