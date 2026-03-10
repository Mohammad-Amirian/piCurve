#' Estimate total column ozone climatologically
#'
#' Estimates total column ozone from latitude, longitude, and date using the
#' climatological expression of Van Heuklon (1979).
#'
#' @param latitude Numeric scalar. Latitude in decimal degrees.
#' @param longitude Numeric scalar. Longitude in decimal degrees.
#' @param date A date coercible to `Date`.
#'
#' @return Numeric scalar giving estimated total column ozone in Dobson Units.
#'
#' @details
#' This function implements a climatological approximation to total column ozone.
#' It can be used when observed ozone depth is unavailable. The returned value is
#' in Dobson Units (DU). If another unit is needed by downstream calculations,
#' conversion should be applied explicitly.
#'
#' @references
#' Van Heuklon, T. K. (1979). Estimating atmospheric ozone for solar radiation
#' models.
#'
#' @examples
#' ozone_depth_piCurve(latitude = 44.65, longitude = -63.57, date = "2024-07-24")
#'
#' @export
ozone_depth_piCurve <- function(latitude, longitude, date) {
    if (!is.numeric(latitude) || length(latitude) != 1L || is.na(latitude)) {
        stop("`latitude` must be a single non-missing numeric value.")
    }
    if (!is.numeric(longitude) || length(longitude) != 1L || is.na(longitude)) {
        stop("`longitude` must be a single non-missing numeric value.")
    }

    date <- as.Date(date)
    if (is.na(date)) {
        stop("`date` must be coercible to class Date.")
    }

    radians <- function(x) x * pi / 180
    E <- as.POSIXlt(date)$yday + 1

    J <- 235
    D <- 0.9865
    G <- 20

    if (latitude > 0) {
        # Northern Hemisphere
        A <- 150
        C <- 40
        FF <- -30
        H <- 3
        beta <- 1.28
        I <- ifelse(longitude > 0, 20, 0)
    } else {
        # Southern Hemisphere
        A <- 100
        C <- 30
        FF <- -152.625
        H <- 20
        beta <- 1.5
        I <- -75
    }

    J + (
        A +
            C * sin(radians(D * (E + FF))) +
            G * sin(radians(H * (longitude + I)))
    ) * (sin(radians(beta * latitude)))^2
}
