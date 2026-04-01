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

#' Estimate UTC hour from local clock time and longitude
#'
#' Converts local clock time to UTC using longitude-based time zone estimation.
#'
#' @param date Date or POSIXct/POSIXt.
#' @param time Character (HH:MM:SS) or POSIXct representing local clock time.
#' @param longitude Longitude in degrees (east positive, west negative).
#'
#' @return Numeric UTC hour (fractional).
#' @keywords internal
.estimate_hour_utc <- function(date, time, longitude) {

    # Parse time
    if (inherits(time, "POSIXt")) {
        hour_local <- lubridate::hour(time) +
            lubridate::minute(time) / 60 +
            lubridate::second(time) / 3600
    } else if (is.character(time)) {
        time_parsed <- lubridate::hms(time)
        hour_local <- lubridate::hour(time_parsed) +
            lubridate::minute(time_parsed) / 60 +
            lubridate::second(time_parsed) / 3600
    } else {
        stop("`time` must be character (HH:MM:SS) or POSIXct/POSIXt.")
    }

    # Theoretical UTC offset (integer timezone approximation)
    utc_offset <- round(longitude / 15)

    hour_utc <- hour_local - utc_offset

    return(hour_utc)
}

#' Local solar time from UTC
#'
#' Converts UTC hour and longitude into local apparent solar time (hours),
#' using an approximate equation of time correction.
#'
#' @param date Date or POSIXct/POSIXt.
#' @param hour_utc Numeric hours since 00:00 UTC.
#' @param longitude Longitude in degrees (east positive, west negative).
#'
#' @return Numeric local solar time in hours.
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
#'
#' @return Numeric hour angle in degrees.
#' @keywords internal
.solar_hour_angle <- function(local_solar_time) {
    15 * (local_solar_time - 12)
}

#' Solar zenith angle
#'
#' Computes solar zenith angle (SZA) in degrees using latitude, longitude,
#' date, and **local clock time**. The function internally estimates UTC time
#' from longitude, so users do not need to perform time zone conversions.
#'
#' @param date Date or POSIXct/POSIXt.
#' @param time Local clock time. Either:
#'   \itemize{
#'     \item Character string "HH:MM:SS"
#'     \item POSIXct/POSIXt object
#'   }
#' @param latitude Latitude in degrees (north positive).
#' @param longitude Longitude in degrees (east positive, west negative).
#'
#' @return Numeric vector of solar zenith angle in degrees in \eqn{[0, 180]}.
#'
#' @details
#' The UTC time is approximated from longitude using:
#' \deqn{UTC\ offset \approx \mathrm{round}(longitude / 15)}
#' This assumes standard time zones and ignores political boundaries and daylight saving time.
#'
#' @examples
#' # Halifax, Canada (~UTC-4), local time 15:00
#' solar_zenith_angle_piCurve(
#'   date = "2026-04-03",
#'   time = "15:00:00",
#'   latitude = 44.65,
#'   longitude = -63.57
#' )
#'
#' @export
solar_zenith_angle_piCurve <- function(date, time, latitude, longitude) {

    date <- as.Date(date)
    # Input checks
    if (!(inherits(date, "Date") || inherits(date, "POSIXct") || inherits(date, "POSIXt"))) {
        stop("`date` must be a Date or POSIXct/POSIXt.")
    }
    if (!is.numeric(latitude) || any(!is.finite(latitude))) {
        stop("`latitude` must be finite numeric.")
    }
    if (!is.numeric(longitude) || any(!is.finite(longitude))) {
        stop("`longitude` must be finite numeric.")
    }

    # Step 1: estimate UTC hour
    hour_utc <- .estimate_hour_utc(date, time, longitude)

    # Step 2: local solar time
    lst <- .local_solar_time_utc(date, hour_utc, longitude)

    # Step 3: hour angle
    h <- .solar_hour_angle(lst)

    # Step 4: declination
    dec <- solar_declination_piCurve(date, hour = hour_utc)

    # Convert to radians
    lat_r <- latitude * pi / 180
    dec_r <- dec * pi / 180
    h_r   <- h * pi / 180

    cos_sza <- sin(lat_r) * sin(dec_r) +
        cos(lat_r) * cos(dec_r) * cos(h_r)

    cos_sza <- pmin(1, pmax(-1, cos_sza))

    acos(cos_sza) * 180 / pi
}
