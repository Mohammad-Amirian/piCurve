#' Estimate sea-surface photosynthetically available radiation (PAR)
#'
#' Estimates clear-sky sea-surface photosynthetically available radiation (PAR)
#' over a user-specified wavelength interval using
#' [sea_surface_irradiance_piCurve()].
#'
#' @param zenith_angle Numeric scalar. Solar zenith angle in degrees
#'   (`0 <= zenith_angle < 90`). Users may estimate solar zenith angle using
#'   [solar_zenith_angle_piCurve()].
#' @param database_info Data frame containing wavelength-specific optical inputs.
#'   Must contain the columns required by [sea_surface_irradiance_piCurve()].
#'   For the built-in dataset, see [GreggCarder1990].
#' @param sample_date A date coercible to `Date`. Used to compute orbital
#'   correction of extraterrestrial irradiance.
#' @param H_oz Numeric scalar. Total column ozone depth in Dobson Units (DU).
#'   Default is `266`.
#' @param data_climato_solar Optional list or one-row data frame of atmospheric
#'   inputs. If `NULL`, default climatological values are used.
#' @param lambda_min Numeric scalar. Lower wavelength bound for PAR integration
#'   in nm. Default is `400`.
#' @param lambda_max Numeric scalar. Upper wavelength bound for PAR integration
#'   in nm. Default is `700`.
#' @param units Character string specifying output units. Must be one of
#'   `"umol photons m^-2 s^-1"` or `"W m^-2"`(default).
#'
#' @return A one-row data frame with columns:
#' \describe{
#'   \item{lambda_min}{Lower wavelength bound of integration (nm).}
#'   \item{lambda_max}{Upper wavelength bound of integration (nm).}
#'   \item{units}{Output PAR units.}
#'   \item{PAR}{Integrated sea-surface PAR in the requested units.}
#' }
#'
#' @details
#' This function first estimates spectral downwelling irradiance at the sea
#' surface using [sea_surface_irradiance_piCurve()], then integrates over the
#' photosynthetically active waveband using the trapezoidal rule.
#'
#' By default, PAR is computed over 400--700 nm.
#'
#' If `units = "W m^-2"`, spectral irradiance is integrated directly:
#' \deqn{
#' PAR = \int_{\lambda_{min}}^{\lambda_{max}} E_d(\lambda)\ d\lambda
#' }
#'
#' If `units = "umol photons m^-2 s^-1"`, spectral irradiance is first converted
#' to photon flux density at each wavelength:
#' \deqn{
#' Q(\lambda) = E_d(\lambda)\frac{\lambda \times 10^{-9}}{h c}
#' \times \frac{10^6}{N_A}
#' }
#' where \eqn{E_d(\lambda)} is spectral irradiance, \eqn{h} is Planck's
#' constant, \eqn{c} is the speed of light, and \eqn{N_A} is Avogadro's number.
#' The converted spectrum is then integrated over wavelength.
#'
#' Because photon energy depends on wavelength, conversion to quantum units must
#' be performed spectrally before integration rather than by applying a single
#' factor to integrated energy-based PAR.
#'
#' The dataset [GreggCarder1990] provides wavelength-specific coefficients
#' consistent with the Gregg and Carder (1990) formulation.
#'
#' @source
#' Gregg, W. W., & Carder, K. L. (1990). A simple spectral solar irradiance
#' model for cloudless maritime atmospheres. \emph{Journal of Geophysical Research}.
#'
#' @examples
#' # Example using Gregg & Carder (1990) dataset
#' PAR_sea_surface_piCurve(
#'   zenith_angle = 30,
#'   database_info = GreggCarder1990,
#'   sample_date = "2024-07-24"
#' )
#'
#' # Return PAR in energy units
#' PAR_sea_surface_piCurve(
#'   zenith_angle = 30,
#'   database_info = GreggCarder1990,
#'   sample_date = "2024-07-24",
#'   units = "umol photons m^-2 s^-1"
#' )
#'
#' # Restrict integration bounds
#' PAR_sea_surface_piCurve(
#'   zenith_angle = 30,
#'   database_info = GreggCarder1990,
#'   sample_date = "2024-07-24",
#'   lambda_min = 450,
#'   lambda_max = 650
#' )
#'
#' # Example based on coordinates and time
#' DATE <- as.Date("2024-03-20")
#' zn <- solar_zenith_angle_piCurve(
#'   DATE,
#'   hour_utc = 12,
#'   latitude = 0,
#'   longitude = 0
#' )
#'
#' PAR_sea_surface_piCurve(
#'   zenith_angle = zn,
#'   database_info = GreggCarder1990,
#'   sample_date = DATE
#' )
#'
#' @export
PAR_sea_surface_piCurve <- function(
        zenith_angle,
        database_info,
        sample_date,
        H_oz = 266,
        data_climato_solar = NULL,
        lambda_min = 400,
        lambda_max = 700,
        units = c("W m^-2", "umol photons m^-2 s^-1")
) {

    units <- match.arg(units)

    if (!is.numeric(lambda_min) || length(lambda_min) != 1L || is.na(lambda_min)) {
        stop("`lambda_min` must be a single non-missing numeric value.")
    }
    if (!is.numeric(lambda_max) || length(lambda_max) != 1L || is.na(lambda_max)) {
        stop("`lambda_max` must be a single non-missing numeric value.")
    }
    if (lambda_min >= lambda_max) {
        stop("`lambda_min` must be smaller than `lambda_max`.")
    }

    spectral_out <- sea_surface_irradiance_piCurve(
        zenith_angle = zenith_angle,
        database_info = database_info,
        sample_date = sample_date,
        H_oz = H_oz,
        data_climato_solar = data_climato_solar
    )

    spectral_out <-
        spectral_out[
            spectral_out$wavelength >= lambda_min &
                spectral_out$wavelength <= lambda_max,
            ,
            drop = FALSE
            ]

    if (nrow(spectral_out) < 2) {
        stop("Not enough wavelength points within the requested integration range.")
    }

    spectral_out <- spectral_out[order(spectral_out$wavelength), , drop = FALSE]

    x <- spectral_out$wavelength
    y <- spectral_out$irradiance

    if (units == "W m^-2") {
        PAR <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    } else {
        h <- 6.62607015e-34
        c <- 2.99792458e8
        N_A <- 6.02214076e23

        lambda_m <- x * 1e-9
        q_lambda <- y * lambda_m / (h * c) * 1e6 / N_A
        PAR <- sum(diff(x) * (head(q_lambda, -1) + tail(q_lambda, -1)) / 2)
    }

    data.frame(
        lambda_min = lambda_min,
        lambda_max = lambda_max,
        units = units,
        PAR = PAR
    )
}

