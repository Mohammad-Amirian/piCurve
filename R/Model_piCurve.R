
#' @title Formulate PI Curve
#' @description
#' Evaluate one of a large number of photosynthesis-irradiance (PI) models as a function of irradiance.
#' Well-established
#' functions such as bi-linear, exponential, hyperbolic tangent, second order polynomial,
#' generalized rectangular hyperbola, and their appropriate combinations are provided.
#'
#' @param parameters Vector -- containing the values listed below:
#' \itemize{
#'      \item{\code{Pmax} \eqn{\hspace{0.1cm}}: }{Maximum photosynthetic rate  (sign +),}
#'      \item{\code{alpha}: }{Slope of the light-saturation curve at low light (sign +),}
#'      \item{\code{beta}  \eqn{\hspace{0.15cm}}: }{Photoinhibition parameter (sign +),}
#'      \item{\code{R}  \eqn{\hspace{0.7cm}}: }{Dark reaction rate (sign + OR -),}
#'      \item{\code{shape}  \eqn{\hspace{0cm}}: }{Shape parameter (sign +) -- Only required for some of the models. See \samp{Details}}.
#' }
#' @param model_name String -- which model? (List of available models is given in \samp{details})
#' @param data Vector -- containing the irradiance profile.
#'
#' @return A vector of predicted photosynthetic rate at each irradiance using the selected function.
#'
#' @export
#' @references
#' Amirian M.M., V Finkel Z., Devred E., Irwin A.J.,
#' “\emph{A new parameterization of photoinhibition for phytoplankton},
#' arXiv (2024) 1–33. 10.48550/arXiv.2412.17923.
#'
#' @details
#' This package comprises three distinct model frameworks for formulating the PI curve:
#' light-limited, light-saturating, and photoinhibition models. Each framework is detailed below
#'
#' ==============================
#'
#' \strong{Table 1}. Model used to formulate light-limited PI curve.
#' | \strong{Name} \eqn{\hspace{0.5cm}} | \strong{Model} | \strong{Function Type} | \strong{References} |
#' | --- | --- |  --- |  --- |
#' | lm | \eqn{P^B = \alpha I  \hspace{1cm}} | Linear | \emph{Linear Regression} |
#'
#' ==============================
#'
#' \strong{Table 2}. Models used to formulate light-saturating PI curve.
#' | \strong{Name} \eqn{\hspace{0.5cm}} | \strong{Model} | \strong{Function Name} | \strong{References} |
#' | --- | --- |  --- |  --- |
#' | LS1 | \eqn{P^B = P^B_{\max}\dfrac{I+I_\alpha - |I-I_\alpha|}{2I_\alpha} \hspace{1cm}} | Blackman-pw | \emph{Blackman 1905} |
#' | LS2 | \eqn{P^B = P^B_{\max}\dfrac{I}{I + I_\alpha}} | Baly-RH | \emph{Baly 1935} |
#' | LS3 | \eqn{P^B = P^B_{\max}\dfrac{I}{\sqrt{I^2 + I^2_\alpha}}} | Smith-RH | \emph{Smith 1936} |
#' | LS4 | \eqn{P^B = P^B_{\max}(1-e^{-I/I_\alpha})} | Webb-exp | \emph{Webb et al. 1974} |
#' | LS5 | \eqn{P^B = P^B_{\max} \tanh \left( \dfrac{I}{I_\alpha} \right)} | Jassby-tanh |  \emph{Jassby et al. 1976} |
#' | LS6 | \eqn{P^B = \dfrac{P^B_{\max}}{2 \theta} \left[\mathcal{I} -\sqrt{\mathcal{I}^2 - 4 \theta(\mathcal{I} - 1)} \right]} \eqn{\hspace{0.5cm}} | Prioul-nonRH | \emph{Prioul et al. 1977} |
#' | LS7 | \eqn{P^B = P^B_{\max} \dfrac{I}{(I^b+I^b_\alpha)^{1/b}}} | Bannister-gRH \eqn{\hspace{0.5cm}} | \emph{Bannister 1979} |
#'
#' ==============================
#'
#' \strong{Table 3}. Models used to formulate photoinhibition PI curve.
#' | \strong{Name} \eqn{\hspace{0.5cm}} | \strong{Model} | \strong{Function Name} | \strong{References} |
#' | --- | --- |  --- |  --- |
#' | Ph01 | \eqn{P^B = P^B_s~ (I/I_{\alpha}^s) ~e^{1-( I/ I_{\alpha}^s)} \hspace{1cm}} | Steele-exp | \emph{Steele 1962} |
#' | Ph02 | \eqn{P^B = P^B_s~ \dfrac{(I/I_{\alpha}^s)}{ ~ I^2 / (I_{\alpha}^s I_{\beta}^s) + (I/I_{\alpha}^s) + 1}} | Peeters-rational | \emph{Peeters et al. 1978} |
#' | Ph03 | \eqn{P^B = P^B_s~ (1-e^{-I/I_{\alpha}^s})e^{-I/I_{\beta}^s} \hspace{1cm}} | Platt-exp | \emph{Platt et al. 1980} |
#' | Ph04 | \eqn{P^B = P^B_s~ \tanh \left(I/I_{\alpha}^s\right) ~ e^{-I/I_{\beta}^s} \hspace{1cm}} | Neale-exp | \emph{Neale et al. 1986} |
#' | Ph05 | \eqn{P^B = P^B_s~ \dfrac{I}{I + I_{\alpha}^s}  ~ e^{-I/I_{\beta}^s} \hspace{1cm}} | Baly-exp | \emph{...} |
#' | Ph06 | \eqn{P^B = P^B_s~ \dfrac{I}{\sqrt{I^2 + (I_{\alpha}^s)^2}}  ~ e^{-I/I_{\beta}^s} \hspace{1cm}} | Smith-exp | \emph{...} |
#' | Ph07 | \eqn{P^B = P^B_s~ \dfrac{I}{(I^b + (I_{\alpha}^s)^b)^{1/b}} ~ e^{-I/I_{\beta}^s} \hspace{1cm}} | Bannister-exp | \emph{...} |
#' | Ph08 | \eqn{P^B = P^B_s~ \dfrac{P^B_{s}}{2 \theta} \left[ \mathcal{I} - \sqrt{\mathcal{I}^2 - 4 \theta(\mathcal{I} - 1)} \right] ~ e^{-I/I_{\beta}^s} \hspace{1cm}} | Prioul-exp | \emph{...} |
#' | Ph09 | \eqn{P^B = \begin{cases} \alpha I & \hspace{0.5cm} I \leq P^B_{\max} / \alpha \\ & \\ P^B_{\max} & \hspace{0.5cm} P^B_{\max} / \alpha < I \leq P^B_{\max} / \beta\\ & \\ -\beta I & \hspace{0.5cm} I > P^B_{\max} / \beta \end{cases} \hspace{1cm}} | Extended-Blackman | \emph{...} |
#' | Ph10 | \eqn{P^B = P^B_{\max} \tanh \left( \dfrac{I}{I_{\alpha}} \right) \tanh{ \left[ \left( \dfrac{I_{\beta}}{I} \right)^{\cosh^2(1)} \right] }  \hspace{1cm}} | Double-tanh | \emph{Amirian et al. 2024} |
#' | Ph11 | \eqn{P^B = P^B_{\max} \tanh \left( \dfrac{I}{I_{\alpha}} \right) \tanh{ \left[ \left( \dfrac{I_{\beta}}{I} \right)^{\gamma} \right] }  \hspace{1cm}} | Double-tanh-shp | \emph{Amirian et al. 2024} |
#' | Ph12 | \eqn{P^B = P^B_{\max}  \left[1 - \exp{\left( -\dfrac{I}{I_{\alpha}} \right)} \right] \tanh{ \left[ \left( \dfrac{I_{\beta}}{I} \right)^{\cosh^2(1)} \right] }  \hspace{1cm}} | Exp-tanh | \emph{Amirian et al. 2024} |
#' | Ph13 | \eqn{P^B = P^B_{\max} \left[1 - \exp{\left( -\dfrac{I}{I_{\alpha}} \right)} \right] \tanh{ \left[ \left( \dfrac{I_{\beta}}{I} \right)^{\gamma} \right] }  \hspace{1cm}} | Exp-tanh-shp | \emph{Amirian et al. 2024} |
#' | Ph14 | \eqn{P^B = P^B_{\max}  \tanh{\left(\dfrac{I}{I_{\alpha}}\right)} \left[1 - \exp{\left(-\dfrac{I_{\beta}}{I}\right)} \right]    \hspace{1cm}} | Tanh-rcp_exp | \emph{Amirian et al. 2024} |
#' | Ph15 | \eqn{P^B = P^B_{\max} \left[1 - \exp{\left(-\dfrac{I}{I_{\alpha}}\right)} \right] \left[1 - \exp{\left(-\dfrac{I_{\beta}}{I}\right)} \right]  \hspace{1cm}} | Exp_rcp_exp | \emph{Amirian et al. 2024} |
#' | Ph16 | \eqn{P^B = \dfrac{P^B_{\max}}{2 \theta} \left[1 + \theta_\beta ~ \mathcal{I}_\alpha -\sqrt{ \theta_\beta ~  \mathcal{I}_\alpha^2 - 4 \theta ~ \mathcal{I}_\alpha + 1} \right]} \eqn{\hspace{0.5cm}}  | Fasham-nonRH | \emph{Fasham et al. 1983}
#'
#' In equations LS6 and Ph08, \eqn{b} is a shape parameter, \eqn{\theta = \left[ 1 + \exp{(-b)} \right]^{-1} } is
#' a sigmoid function setting \eqn{0< \theta <1},
#' \eqn{\theta_\beta = \theta + (1-\theta)\exp{(\beta I)}},
#' \eqn{\mathcal{I} = \left({I}/{I_\alpha} + 1\right) },
#' and \eqn{\mathcal{I}_\alpha = I / I_\alpha}.
#' In the above tables also, \eqn{I_{\alpha} = P_{\max}^B/ \alpha}, \eqn{I^s_{\alpha} = P_{s}^B/ \alpha}, \eqn{I_{\beta} = P_{\max}^B/ \beta}, and \eqn{I^s_{\beta} = P_{s}^B/ \beta}.
#'
#' ==============================
#'
#' Note that the default model_name for light-saturating and photoinhibition data types is
#' Jassby-tanh (LS5) and double-tanh (Ph10) model, resepctively.
#' @import tibble

#' @examples
#' # model parameters
#' params <- c(Pmax = 20, alpha = 0.6, R = 0)
#'
#' # generate an irradiance profile
#' df <- tibble::tibble(I = seq(0, 300, length = 25))
#'
#' # compute the photosynthetic rate using Jassby-tanh (LS5) model
#' Model_piCurve(parameters = params, model_name = "tanh", data = df)
#'
#' # compute the photosynthetic rate using double-tanh (Ph10) model
#' Model_piCurve(parameters = c(params, beta = 0.3), model_name = "2tanh", data = df)
#'
Model_piCurve <-
    function(parameters,
             model_name,
             data
             ){

        # data format check ----
        if (!all("I" %in% names(data))) {
            stop(print("Input data frame must contain 'I' columns!"))
        }

        # data type = ? ----
        data_type <-
            ifelse(
                nameFormat(model_name) %in% Pool_eqName_lm, "ll", # light-limited
                ifelse(nameFormat(model_name) %in% Pool_eqName_ls,
                       "ls", # light-saturating
                       "ph" # photoinhibition
                       )
                )



        # select the model specified by the user
        equation <- Model_setup(which_model = model_name, data_type)

        # compute the primary production profile and return it
        return(equation(parameters = parameters, data))

        }




