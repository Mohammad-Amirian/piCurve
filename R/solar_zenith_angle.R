#' Equation of time (minutes)
#'
#' Approximate equation of time (EoT) in minutes as a function of day-of-year.
#'
#' @param date Date or POSIXct/POSIXt.
#' @return Numeric vector of EoT in minutes.
#' @keywords internal
.equation_of_time <- function(date) {
    doy <- lubridate::yday(date)
    B <- (2 * pi / 365) * (doy - 81)
    9.87 * sin(2 * B) - 7.53 * cos(B) - 1.5 * sin(B)
}

#' Local solar time from UTC
#'
#' Converts UTC hour and longitude into local apparent solar time (hours),
#' using an approximate equation of time correction.
#'
#' @param date Date or POSIXct/POSIXt (used for day-of-year in EoT).
#' @param hour_utc Numeric hours since 00:00 UTC (may be fractional).
#' @param longitude Longitude in degrees (east positive, west negative).
#' @return Numeric local solar time in hours (not wrapped).
#' @keywords internal
.local_solar_time_utc <- function(date, hour_utc, longitude) {
    eot_min <- .equation_of_time(date)
    hour_utc + (longitude / 15) + (eot_min / 60)
}

#' Solar hour angle
#'
#' Hour angle in degrees. Negative in the morning, positive in the afternoon.
#'
#' @param local_solar_time Local solar time in hours.
#' @return Numeric hour angle in degrees.
#' @keywords internal
.solar_hour_angle <- function(local_solar_time) {
    15 * (local_solar_time - 12)
}

#' Solar zenith angle
#'
#' Computes solar zenith angle (SZA) in degrees using latitude, longitude, date,
#' and UTC hour, with an approximate equation of time correction.
#'
#' @param date Date or POSIXct/POSIXt.
#' @param hour_utc Numeric hours since 00:00 UTC (may be fractional). Length 1
#'   or same length as `date`.
#' @param latitude Latitude in degrees (north positive).
#' @param longitude Longitude in degrees (east positive, west negative).
#'
#' @return Numeric vector of solar zenith angle in degrees in \eqn{[0, 180]}.
#'
#' @examples
#' # Equator at equinox: near zenith at local noon
#' solar_zenith_angle_piCurve(as.Date("2024-03-20"), hour_utc = 12, latitude = 0, longitude = 0)
#'
#' @export
solar_zenith_angle_piCurve <- function(date, hour_utc, latitude, longitude) {

    # Basic checks (keep these lightweight if you prefer)
    if (!(inherits(date, "Date") || inherits(date, "POSIXct") || inherits(date, "POSIXt"))) {
        stop("`date` must be a Date or POSIXct/POSIXt.")
    }
    if (!is.numeric(hour_utc) || any(!is.finite(hour_utc))) stop("`hour_utc` must be finite numeric.")
    if (!is.numeric(latitude) || any(!is.finite(latitude))) stop("`latitude` must be finite numeric.")
    if (!is.numeric(longitude) || any(!is.finite(longitude))) stop("`longitude` must be finite numeric.")

    # Local solar time (hours)
    lst <- .local_solar_time_utc(date, hour_utc, longitude)

    # Hour angle (deg)
    h <- .solar_hour_angle(lst)

    # Declination (deg)
    dec <- solar_declination_piCurve(date, hour = hour_utc)

    # Radians
    lat_r <- latitude * pi / 180
    dec_r <- dec * pi / 180
    h_r   <- h * pi / 180

    cos_sza <- sin(lat_r) * sin(dec_r) + cos(lat_r) * cos(dec_r) * cos(h_r)

    # Clamp for numerical safety
    cos_sza <- pmin(1, pmax(-1, cos_sza))

    acos(cos_sza) * 180 / pi
}
