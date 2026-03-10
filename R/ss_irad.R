#' Estimate clear-sky sea-surface spectral irradiance
#'
#' Estimates downwelling clear-sky spectral irradiance at the sea surface using
#' a Gregg and Carder (1990) parameterization for cloudless maritime atmospheres.
#'
#' @param zenith_angle Numeric scalar. Solar zenith angle in degrees
#'   (`0 <= zenith_angle < 90`). Users may estimate solar zenith angle using
#'   [solar_zenith_angle_piCurve()].
#' @param database_info Data frame containing wavelength-specific optical inputs.
#'   Must contain the following columns:
#'   \describe{
#'     \item{lambda}{Wavelength in vacuo (nm).}
#'     \item{H_o}{Spectral extraterrestrial solar irradiance (W m\eqn{^{-2}} nm\eqn{^{-1}}).}
#'     \item{a_oz}{Ozone absorption coefficient (cm\eqn{^{-1}}).}
#'     \item{a_w}{Water vapor absorption coefficient (cm\eqn{^{-1}}).}
#'     \item{a_o}{Oxygen absorption coefficient (cm\eqn{^{-1}}).}
#'   }
#' @param sample_date A date coercible to `Date`. Used to compute orbital
#'   correction of extraterrestrial irradiance.
#' @param H_oz Numeric scalar. Total column ozone depth in Dobson Units (DU).
#'   Default is `266`. If not supplied, users may estimate ozone depth using
#'   [ozone_depth_piCurve()].
#' @param data_climato_solar Optional list or one-row data frame of atmospheric
#'   inputs (e.g., pressure, water vapor, aerosol properties). If `NULL`,
#'   default climatological values are used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{wavelength}{Wavelength (nm).}
#'   \item{irradiance}{Estimated downwelling spectral irradiance at the sea
#'   surface (W m\eqn{^{-2}} nm\eqn{^{-1}}).}
#' }
#'
#' @details
#' This function implements a spectral clear-sky irradiance model following
#' Gregg and Carder (1990), combining atmospheric transmittance due to:
#'
#' \itemize{
#'   \item Rayleigh scattering
#'   \item Oxygen absorption
#'   \item Water vapor absorption
#'   \item Ozone absorption
#'   \item Aerosol scattering and absorption
#'   \item Air–sea surface reflectance
#' }
#'
#' Spectral irradiance is computed as a function of wavelength and solar
#' zenith angle, with extraterrestrial irradiance corrected for Earth–Sun
#' distance based on day of year.
#'
#' \strong{NOTE: Unit consistency is required:}
#' \itemize{
#'   \item Wavelength (`lambda`) must be provided in nanometers (nm).
#'   \item Extraterrestrial irradiance (`H_o`) must be in W m\eqn{^{-2}} nm\eqn{^{-1}}.
#'   \item Absorption coefficients (`a_oz`, `a_w`, `a_o`) must be in cm\eqn{^{-1}}.
#'   \item Ozone depth (`H_oz`) must be in Dobson Units (DU).
#' }
#'
#' The function returns spectral irradiance. Photosynthetically available
#' radiation (PAR) can be obtained by integrating output irradiance over
#' 400–700 nm.
#'
#' The dataset [GreggCarder1990] provides a standard set of coefficients
#' consistent with this formulation.
#'
#' @source
#' Gregg, W. W., & Carder, K. L. (1990). A simple spectral solar irradiance
#' model for cloudless maritime atmospheres. \emph{Journal of Geophysical Research}.
#'
#' @examples
#' # Example with synthetic coefficients
#' db <- data.frame(
#'   lambda = c(400, 500, 600),
#'   a_o = c(0.01, 0.02, 0.01),
#'   a_w = c(0.00, 0.01, 0.02),
#'   a_oz = c(0.03, 0.02, 0.01),
#'   H_o = c(1.8, 1.9, 1.7)
#' )
#'
#' sea_surface_irradiance_piCurve(
#'   zenith_angle = 30,
#'   database_info = db,
#'   sample_date = "2024-07-24"
#' )
#'
#' # Example using Gregg & Carder (1990) dataset
#' sea_surface_irradiance_piCurve(
#'   zenith_angle = 30,
#'   database_info = GreggCarder1990,
#'   sample_date = "2024-07-24"
#' )
#' @export
sea_surface_irradiance_piCurve <- function(
        zenith_angle,
        database_info,
        sample_date,
        H_oz = 266,
        data_climato_solar = NULL
) {
    req_cols <- c("lambda", "a_o", "a_w", "a_oz", "H_o")

    if (!is.numeric(zenith_angle) || length(zenith_angle) != 1L || is.na(zenith_angle)) {
        stop("`zenith_angle` must be a single non-missing numeric value.")
    }
    if (zenith_angle < 0 || zenith_angle >= 90) {
        stop("`zenith_angle` must be in [0, 90).")
    }
    if (!is.data.frame(database_info)) {
        stop("`database_info` must be a data frame.")
    }
    if (!all(req_cols %in% names(database_info))) {
        stop("`database_info` must contain columns: ",
             paste(req_cols, collapse = ", "), ".")
    }

    sample_date <- as.Date(sample_date)
    if (is.na(sample_date)) {
        stop("`sample_date` must be coercible to class Date.")
    }

    if (is.null(data_climato_solar)) {
        P <- 1010.5
        WV <- 3.25
        beta <- 0.27
        alpha <- 0.5
        AM <- 5.5
        RH <- 68.5
        W <- 2.83
    } else {
        P <- data_climato_solar$P
        WV <- data_climato_solar$WV
        alpha <- data_climato_solar$alpha
        beta <- (3.912 / data_climato_solar$V) * (database_info$lambda)^alpha
        AM <- data_climato_solar$AM
        RH <- data_climato_solar$RH
        W <- data_climato_solar$W
        H_oz <- data_climato_solar$O3
    }

    mass_val <- atmospheric_mass_piCurve(zenith_angle, non_stnd_pressure = P)
    m <- unname(mass_val["M_theta"])
    m_prime <- unname(mass_val["M_theta_prime"])
    m_ozone <- unname(mass_val["M_ozone"])

    T_r <- trans_rayleigh(corrected_air_mass = m_prime, wavelength = database_info$lambda)

    T_O2_W <- trans_O2_WaterVapor(
        O2abs_coef = database_info$a_o,
        WVabs_coef = database_info$a_w,
        air_mass = m,
        corrected_air_mass = m_prime,
        WV = WV
    )

    T_o <- T_O2_W$trans_O2
    T_w <- T_O2_W$trans_wv

    T_oz <- trans_ozone(
        ozone_abs_coef = database_info$a_oz,
        ozone_mass = m_ozone,
        ozone_height = H_oz
    )

    aerosol_out <- aerosol_transmittance(
        wavelength = database_info$lambda,
        zenith_angle = zenith_angle,
        air_mass = m,
        beta = beta,
        alpha = alpha,
        AM = AM,
        RH = RH
    )

    T_a <- aerosol_out$trans_aerosol
    T_aa <- aerosol_out$Taa
    T_as <- aerosol_out$Tas

    rho <- surface_reflectance(zenith_angle = zenith_angle, W = W)
    rho_d <- unname(rho["rho_d"])
    rho_s <- unname(rho["rho_s"])

    Fo <- extraterrestrial_irradiance_corrected(date = sample_date, mean_ext_radiance = database_info$H_o)

    incoming_term <- Fo * cos(zenith_angle * pi / 180) * T_o * T_w * T_oz
    Ia <- incoming_term * T_aa * (T_r^1.5) * (1 - T_as)
    Ir <- incoming_term * T_aa * (1 - T_r^0.95)

    Edd <- incoming_term * T_r * T_a * (1 - rho_d)
    Eds <- (Ia + Ir) * (1 - rho_s)
    Ed <- Edd + Eds

    data.frame(
        wavelength = database_info$lambda,
        irradiance = Ed
    )
}
