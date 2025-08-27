#' Photosynthesis–Irradiance Dataset
#'
#' A dataset containing 8 independent PI incubation samples for model testing and validation.
#'
#' @format A data frame with `n` rows and 4 variables:
#' \describe{
#'   \item{I}{Irradiance (W m⁻²)}
#'   \item{P}{Photosynthesis rate (mol C (mg chl-a)⁻¹ h⁻¹)}
#'   \item{pi_number}{Unique ID for each PI incubation}
#'   \item{data_type}{Type of data: "light-saturated" or "photoinhibition"}
#' }
#'
#' @source Amirian et al. (2025); raw data in `"data-raw/piDataSet.csv"`
#'
#' @references
#' M. Amirian, M., Finkel, Z.V., Devred, E., Irwin, A.J.
#' "\emph{Parameterization of photoinhibition for phytoplankton}".
#' Commun Earth Environ 6, 707 (2025). https://doi.org/10.1038/s43247-025-02686-3
#'
#' @details
#' For inquiries or issues, please include the corresponding `pi_number` and contact \verb{m.amirianmalotb@dal.ca}.

#' @examples
#' head(piDataSet, 3)
"piDataSet"
