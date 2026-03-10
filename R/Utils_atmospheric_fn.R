# Utility functions for atmospheric
#' @keywords internal
#' @noRd
#'
fresnel_reflectance <- function(zenith_angle) {
    if (!is.numeric(zenith_angle) || length(zenith_angle) != 1L || is.na(zenith_angle)) {
        stop("`zenith_angle` must be a single non-missing numeric value.")
    }
    if (zenith_angle < 0 || zenith_angle >= 90) {
        stop("`zenith_angle` must be in [0, 90).")
    }

    n_air <- 1.0
    n_water <- 1.33
    n <- n_water / n_air

    if (abs(zenith_angle) < 1e-10) {
        return(((n_air - n_water) / (n_air + n_water))^2)
    }

    theta_i <- zenith_angle * pi / 180
    theta_t <- asin(sin(theta_i) / n)

    R_s <- (sin(theta_i - theta_t) / sin(theta_i + theta_t))^2
    R_p <- (tan(theta_i - theta_t) / tan(theta_i + theta_t))^2

    (R_s + R_p) / 2
}



# calculate atmospheric transmittance due to Rayleigh scattering ----
trans_rayleigh <-
    function(corrected_air_mass, wavelength) {
        if (any(wavelength <= 0, na.rm = TRUE)) {
            stop("`wavelength` must contain positive values.")
        }
        M_theta_prime <- corrected_air_mass
        lambda <- wavelength
        exp(-M_theta_prime / (115.6406 * lambda ^ 4 - 1.335 * lambda ^ 2))
    }

# calculate O2 & water vapor absorption ----
# WV: total precipitate WV (ranges from 1.8 to 4.7 cm) => mean(wv) = 3.25
trans_O2_WaterVapor <-
    function(O2abs_coef,
             WVabs_coef,
             air_mass,
             corrected_air_mass,
             WV  = 3.25){

        a_o <- O2abs_coef
        a_w <- WVabs_coef
        M_theta <- air_mass
        M_theta_prime <- corrected_air_mass

        # transmittance due to O2 absorption
        O2_par <- a_o * M_theta_prime
        trans_O2 <- exp(-1.41 * O2_par / (1 + 118.3 * O2_par)^(0.45))

        # transmittance due to water vapor

        wv_par <- a_w * WV * M_theta
        trans_wv <- exp( (-0.2385 * wv_par) / (1 + 20.07 * wv_par)^(0.45) )

        return(list(trans_O2 = trans_O2, trans_wv = trans_wv))
    }

# calculate ozone absorption ----
trans_ozone <-
    function(ozone_abs_coef, ozone_mass, ozone_height){
        if (any(ozone_height < 0, na.rm = TRUE)) {
            stop("`ozone_height` must be non-negative.")
        }
        a_oz <- ozone_abs_coef
        M_oz <- ozone_mass
        # Convert ozone depth from Dobson Units to atm-cm equivalent scale
        H_oz <- ozone_height * 1e-3

        exp(-a_oz * M_oz * H_oz)
    }


# aerosol optical depth ----
aerosol_transmittance <-
    function(wavelength, zenith_angle, air_mass,
             beta = 0.27, alpha = 0.5, AM = 5.5, RH = 68.5){
        # beta: turbidity_coef
        # alpha: Agn_exponent

        M_theta <- air_mass
        lambda <- wavelength
        theta_radians <- zenith_angle * pi / 180

        # units: turbidity_coef is Km-1 & wavelength is nm

        AOD <- beta * (lambda / 550)^( -alpha )

        Ta <- # trans_aerosol
            exp( -AOD * M_theta)

        # scattering aerosol
        ave_cos <- -0.1417 * alpha + 0.82

        B3 <- log(1-ave_cos)
        B2 <- B3 * (0.0783 + B3 * (-0.3824 - 0.5874 * B3))
        B1 <- B3 * (1.459 + B3 * (0.1595 + 0.4129 * B3))

        Fa <- 1 - 0.5 * exp((B1 + B2 * cos(theta_radians)) * cos(theta_radians))

        # single-scattering albedo
        Wa <- (-0.0032 * AM + 0.972) * exp(3.06 * 1e-04 * RH)


        # transmittace after aerosol
        Taa <- exp( -(1 - Wa) * AOD * M_theta)

        # transmittace after aerosol
        Tas <- exp( -Wa * AOD * M_theta)

        return(list(trans_aerosol = Ta,
                 Fa = Fa,
                 Wa = Wa,
                 Taa = Taa,
                 Tas = Tas))
    }

surface_reflectance <- function(zenith_angle, W = 2.83){

    theta <- zenith_angle

    if (W > 2) {
        if (theta >= 40) {
            b <- -7.14 * 1e-4 * W + 0.0618
            rho_dsp <- 0.0253 * exp(b * (theta - 40))
        } else {
            rho_dsp <- 0
        }

    } else {
        rho_dsp <- fresnel_reflectance(theta)
    }

    rho_a <- 1.2 * 1000 # air density (g m-3)

    if (!is.numeric(W) || length(W) != 1L || is.na(W) || W < 0) {
        stop("`W` must be a single non-missing numeric value >= 0.")
    }

     # drag Coef
    if (W == 0) {
        Cd <- (0.62) * 1e-3
    } else {
        Cd <- ifelse(W > 7,
                     (0.49 + 0.065 * W) * 1e-3,
                     (0.62 + 1.56 * (1 / W)) * 1e-3)
    }

    # other coef values
    D1 <- 2.2 * 1e-5
    D2 <- 4 * 1e-4
    D3 <- 4.5 * 1e-5
    D4 <- D2/10

    rho_f <- ifelse(W <= 4, 0,
                    ifelse(W > 7,
                           (D3 * rho_a * Cd - D4) * W ^ 2,
                           D1 * rho_a * Cd * W^2 - D2))

    rho_ssp <- ifelse(W > 4, 0.057, 0.066)

    rho_d <- rho_dsp + rho_f
    rho_s <- rho_ssp + rho_f

    return(c(rho_d = rho_d, rho_s = rho_s))
}


# extraterrestrial_radiance_corrected
extraterrestrial_irradiance_corrected <-
    function(date, mean_ext_radiance){
        H_o <- mean_ext_radiance

        # Orbital eccentricity
        e_orbital <- 0.0167

        # Extract the day of the year
        D <- lubridate::day(date)
        # Utility functions for atmospheric
        D_adj <- (D - 3) / 365

        orb_eccentricity_term <-
            (1 + e_orbital * cos(2 * pi * D_adj))^2

        return(H_o * orb_eccentricity_term)
    }

