#' Solar declination angle
#'
#' Computes the solar declination angle (Sun's latitude relative to the Earth's
#' equatorial plane) for a given date and fractional hour.
#'
#' @details
#' This function uses a common low-order solar position approximation:
#' it computes the Julian Day (Gregorian calendar), days since J2000.0,
#' the Sun's mean longitude and anomaly, the ecliptic longitude, and then
#' converts to declination using the obliquity of the ecliptic.
#'
#' Returned values are in **degrees**.
#'
#' @param date A `Date` or `POSIXct` vector. If `POSIXct`, the time-of-day in
#'   `date` is **not used** unless you explicitly pass `hour`.
#' @param hour Numeric. Hours since midnight (typically UTC), may be fractional
#'   (e.g., `13.5` for 13:30). Length 1 or `length(date)`.
#'
#' @return A numeric vector of solar declination angles in degrees, same length
#'   as `date`.
#'
#' @examples
#' # Around the March equinox declination is near 0 deg:
#' solar_declination_piCurve(as.Date("2024-03-20"), hour = 3.1)
#'
#' # June solstice is near +23.44 deg:
#' solar_declination_piCurve(as.Date("2024-06-21"))
#'
#' # Vectorized dates:
#' solar_declination_piCurve(as.Date(c("2024-03-20", "2024-06-21", "2024-12-21")))
#'
#' # Include fractional hour:
#' solar_declination_piCurve(as.Date("2024-06-21"), hour = 12)
#'
#' @export
solar_declination_piCurve <- function(date, hour = 0) {

    if (!(inherits(date, "Date") || inherits(date, "POSIXct") || inherits(date, "POSIXt"))) {
        stop("`date` must be a Date or POSIXct/POSIXt.")
    }

    if (!(is.numeric(hour) && all(is.finite(hour)))) {
        stop("`hour` must be a finite numeric vector (can be fractional).")
    }

    n_date <- length(date)
    n_hour <- length(hour)
    if (!(n_hour == 1L || n_hour == n_date)) {
        stop("`hour` must have length 1 or the same length as `date`.")
    }

    # ---- Julian Day (Gregorian calendar) ----
    day   <- lubridate::day(date)
    month <- lubridate::month(date)
    year  <- lubridate::year(date)

    # Vector-safe adjustment for Jan/Feb
    idx <- month <= 2
    year[idx]  <- year[idx] - 1
    month[idx] <- month[idx] + 12

    A <- floor(year / 100)
    B <- 2 - A + floor(A / 4)

    JD <- floor(365.25 * (year + 4716)) +
        floor(30.6001 * (month + 1)) +
        day + B - 1524.5 + (hour / 24)

    # Days since J2000.0
    n <- JD - 2451545.0

    # Mean longitude (deg)
    L <- (280.46 + 0.9856474 * n) %% 360

    # Mean anomaly (deg)
    g <- (357.528 + 0.9856003 * n) %% 360

    # Ecliptic longitude (deg)
    lambda <- (L + 1.915 * sin(g * pi / 180) + 0.020 * sin(2 * g * pi / 180)) %% 360

    # Obliquity (deg)
    epsilon <- 23.439 - 0.0000004 * n

    # Declination (deg)
    delta <- asin(sin(epsilon * pi / 180) * sin(lambda * pi / 180)) * 180 / pi

    delta
}
