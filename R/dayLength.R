#' Calculate astronomical day length
#'
#' Computes day length (photoperiod) for a given date and latitude using
#' solar declination calculated by [solar_declination_piCurve()].
#'
#' @param date A `Date`, `POSIXct`, or `POSIXt` object.
#' @param lat_deg Numeric latitude in degrees. Northern latitudes are positive
#'   and southern latitudes are negative.
#' @param hour Numeric hour of the day used to calculate solar declination.
#'   Default is `12`, representing local noon. Can be fractional.
#'
#' @return Numeric vector of day length in hours.
#'
#' @details
#' Day length is calculated from the sunset hour angle:
#'
#' \deqn{
#' D = \frac{24}{\pi}\cos^{-1}(-\tan(\phi)\tan(\delta))
#' }
#'
#' where \eqn{\phi} is latitude and \eqn{\delta} is solar declination, both in
#' radians. Polar day and polar night cases are handled by constraining the
#' argument of \eqn{\cos^{-1}} to the interval [-1, 1].
#'
#' The `hour` argument is retained for consistency with
#' [solar_declination_piCurve()]. For day-length calculations, `hour = 12` is
#' usually appropriate because within-day changes in solar declination are very
#' small.
#'
#' @examples
#' # Day length at Bedford Basin, Nova Scotia, on summer solstice
#' day_length_piCurve(as.Date("2026-06-21"), lat_deg = 44.69)
#'
#' # Seasonal day length over one year
#' dates <- seq(as.Date("2026-01-01"), as.Date("2026-12-31"), by = "day")
#' dl <- day_length_piCurve(dates, lat_deg = 44.69)
#'
#' @export
day_length_piCurve <- function(date, lat_deg, hour = 12) {

    if (!(inherits(date, "Date") || inherits(date, "POSIXct") ||
          inherits(date, "POSIXt"))) {
        stop("`date` must be a Date or POSIXct/POSIXt.")
    }

    if (!(is.numeric(lat_deg) && all(is.finite(lat_deg)))) {
        stop("`lat_deg` must be a finite numeric vector.")
    }

    if (any(lat_deg < -90 | lat_deg > 90)) {
        stop("`lat_deg` must be between -90 and 90 degrees.")
    }

    if (!(is.numeric(hour) && all(is.finite(hour)))) {
        stop("`hour` must be a finite numeric vector (can be fractional).")
    }

    n_date <- length(date)
    n_lat <- length(lat_deg)
    n_hour <- length(hour)

    if (!(n_lat == 1L || n_lat == n_date)) {
        stop("`lat_deg` must have length 1 or the same length as `date`.")
    }

    if (!(n_hour == 1L || n_hour == n_date)) {
        stop("`hour` must have length 1 or the same length as `date`.")
    }

    delta_deg <- solar_declination_piCurve(date = date, hour = hour)

    phi <- lat_deg * pi / 180
    delta <- delta_deg * pi / 180

    x <- -tan(phi) * tan(delta)

    x <- pmin(pmax(x, -1), 1)

    day_length <- (24 / pi) * acos(x)

    day_length
}
