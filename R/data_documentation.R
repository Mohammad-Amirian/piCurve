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
#' @source Amirian et al. (2024); raw data in `"data-raw/piDataSet.csv"`
#'
#' @references
#' Amirian M.M., Finkel Z.V., Devred E., Irwin A.J. (2024).
#' *A new parameterization of photoinhibition for phytoplankton*. arXiv:2412.17923.
#'
#' @details
#' For inquiries or issues, please include the corresponding `pi_number` and contact \verb{m.amirianmalotb@dal.ca}.

#' @examples
#' head(piDataSet, 3)
"piDataSet"
