
#' @title Formulate PI Curve
#' @description
#' Wrapper to formulate the photosynthesis-irradiance curve using well-established
#' functions such as bi-linear, exponential, hyperbolic tangent, second order polynomial,
#' generalized rectangular hyperbola, and their appropriate combinations.
#'
#' @param parameters Vector -- containing the values listed below:
#' \itemize{
#'      \item{\code{Pmax} \eqn{\hspace{0.1cm}}: }{Maximum photosynthesis rate normalised by Chl_a (sign +),}
#'      \item{\code{alpha}: }{Light-saturation slope at low light level (sign +),}
#'      \item{\code{beta}  \eqn{\hspace{0.15cm}}: }{Photoinhibition rate at high light level (sign +),}
#'      \item{\code{R}  \eqn{\hspace{0.7cm}}: }{Dark reaction parameter (sign + OR -),}
#'      \item{\code{shape}  \eqn{\hspace{0cm}}: }{Shape parameter (sign +) -- Only required for some of the models. See \samp{Details}}.
#' }

#' @param model_name String -- which model? (List of available models is given in \samp{details})
#' @param data Vector -- Containing the irradiance profile.
#'
#' @return A vector that containing the primary production (\verb{PP}) calculated
#' by a given formula provided by user.
#' @export
#' @references
#' Will be add later
#'
#' @details
#' This package includes both commonly-used and recently-developed PI models.
#' The list of models is provided below.
#'
#' \strong{Table 1}. Existing formulas for light-saturation curve (\eqn{I_k = P_{\max}^B/ \alpha}).
#' | \strong{Name} \eqn{\hspace{0.5cm}} | \strong{Model} | \strong{Function Type} | \strong{References} |
#' | --- | --- |  --- |  --- |
#' | Eq1 | \eqn{P^B = \alpha I  \hspace{1cm}} | Linear | \emph{Linear Regression} |
#' | Eq2 | \eqn{P^B = P^B_{\max}\dfrac{I+I_k - |I-I_k|}{2I_k} \hspace{1cm}} | Bilinear | \emph{Blackman 1905} |
#' | Eq3 | \eqn{P^B = P^B_{\max}\dfrac{I}{I + I_k}} | Rectangular Hyperbola | \emph{Baly 1935} |
#' | Eq4 | \eqn{P^B = P^B_{\max}\dfrac{I}{\sqrt{I^2 + I^2_k}}} | Modified Rectangular Hyperbola | \emph{Smith 1936} |
#' | Eq5 | \eqn{P^B = P^B_{\max} \ln \left( \dfrac{2I}{I_k} \right)} | Logarithm | \emph{Talling 1957} |
#' | Eq6 | \eqn{P^B = P^B_{\max}\ln \left( \dfrac{I}{I_k} + \sqrt{1 +  \left(\dfrac{I}{I_k}\right)^2} \right)}| Modified Logarithm | \emph{Vollenweider 1958} |
#' | Eq7 | \eqn{P^B = P^B_{\max}(1-e^{-I/I_k})} | Exponential | \emph{Webb et al. 1974} |
#' | Eq8 | \eqn{P^B = P^B_{\max} \tanh \left( \dfrac{I}{I_k} \right)} | Hyperbolic Tangent|  \emph{Jassby et al. 1976} |
#' | Eq9 | \eqn{P^B = \dfrac{P^B_{\max}}{2 \theta} \left[\mathcal{I} -\sqrt{\mathcal{I}^2 - 4 \theta(\mathcal{I} - 1)} \right]} \eqn{\hspace{0.5cm}} | Non-rectangular hyperbola | \emph{Prioul et al. 1977} |
#' | Eq10 | \eqn{P^B = P^B_{\max} \dfrac{I}{(I^b+I^b_k)^{1/b}}} | Generalized Rectangular Hyperbola \eqn{\hspace{0.5cm}} | \emph{Bannister 1979} |
#'
#' where \eqn{b} is a shape parameter, \eqn{\theta = \dfrac{1}{1 + \exp{(-b)}}} is
#' a sigmoid function setting \eqn{0< \theta <1}, and \eqn{\mathcal{I} = \left(\dfrac{I}{I_k} + 1\right) }.
#'
#' @note
#' The function \samp{Model_piCurve} can handle all combinations of the names
#' listed above, where the input string \verb{model_name}  is concerned.
#' For instance,if any of the following forms are used, the function will set
#'  \verb{model_name = "Eq1_Blackman"}.
#'
#' \itemize{
#'          \item With or without spacing: Eq1Blackman, Eq1 Blackman,
#'          \item Underscore or hyphen spacing: Eq1_Blackman, Eq1-Blackman,
#'          \item Any form of lower/upper case: eq1_blackman, EQ1_blackman,
#'          \item Missing author name: EQ1, Eq1, eQ1, eq1, 1,
#'          \item Missing Eq indicator: blackman, Blackman,
#'          \item Other combination of all the above: Blackmaneq1, BlackmanEq1,
#'          Blackman Eq1, eq1 Blackman, etc.
#' }
#' @import tibble

#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 100, length = 25))
#'
#' # compute the primary production profile using Baly's rectangular hyperbola model
#' Model_piCurve(parameters = params, model_name = "Eq3-Baly", data = df)
#'
#' # compute the primary production profile using Eq10 model
#' Model_piCurve(parameters = c(params, beta = 0.3), model_name = "Eq11_Blackman_extended_1", data = df)
#'
Model_piCurve <-
    function(parameters,
             model_name = c(
                 ## no photoinhibition
                 "Eq1_lm",
                 "Eq2_Blackman",
                 "Eq3_Baly",
                 "Eq4_Smith",
                 "Eq5_Talling",
                 "Eq6_Vollenweider_Ln",
                 "Eq7_Webb",
                 "Eq8_Jassby",
                 "Eq9_Prioul",
                 "Eq10_Bannister",

                 ## with photoinhibition
                 "Eq11_Blackman_extended_1",
                 "Eq12_Blackman_extended_2",
                 "Eq13_Vollenweider",
                 "Eq14_Steele",
                 "Eq15_Steele_modified",
                 "Eq16_Peeters",
                 "Eq17_Platt",
                 "Eq18_Neale",
                 "Eq19_Neale_modified",
                 "Eq20_Ye",
                 "Eq21_Baly_extended",
                 "Eq22_Smith_extended",
                 "Eq23_Prioul_extended",
                 "Eq24_Bannister_extended",
                 "Eq25_Tanh_Tanh_1",
                 "Eq26_Tanh_Tanh_2",
                 "Eq27_Tanh_Tanh_shape",
                 "Eq28_Tanh_Expo",
                 "Eq29_Expo_Tanh_1",
                 "Eq30_Expo_Tanh_2",
                 "Eq31_Expo_Tanh_shape",
                 "Eq32_Expo_Expo",
                 "Eq33_Amirian_Tanh",
                 "Eq34_Amirian_Expo"
                 ),
             data){
        # select the model specified by the user
        equation <- Model_setup(which_model = model_name)

        # compute the primary production profile and return it
        return(equation(parameters = parameters, data))

        }




