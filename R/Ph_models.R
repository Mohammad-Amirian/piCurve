#####################################################################
# Functions used to model pi-curve with with photoinhibition (beta ~= 0)

################## model notations ##################################
# Pmax : maximum photosynthesis rate (sign +)
# alpha: light-saturation slope, low light level (sign +)
# beta : photoinhibition rate at high light level (sign +)
# R    : dark reaction parameter (sign + OR -)
# shape: shape parameter (sign +)

#####################################################################
################## model details ####################################
# Models other than those listed below have 4 parameters: Pmax, alpha, beta, and R.
#### Eq14_Steele has 3 parameters (no beta)

#### Following equations have 5 parameter: the extra param is shape
######## Eq13_Vollenweider
######## Eq23_Prioul_extended
######## Eq24_Bannister_extended
######## Eq27_Tanh_Tanh_shape
######## Eq31_Expo_Tanh_shape

#### Following equations have 5 parameter: the extra param is Ibreak, where photoinhibition kicks
######## Eq12_Blackman_extended_5params
######## Eq18_Neale

#####################################################################
################## formulations #####################################

# Eq11: Blackman model, extended 4 params
Eq11_Blackman_extended_4params <- function(data, parameters){

    # Pmax, alpha, and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <-
        alpha * I + Heaviside(I - Ik) * (Pmax - alpha * I) +
        Heaviside(I - Ikb) * (-beta * (I - Ikb)) + R

    return(PPhat)
}

# Eq12: Blackman model, extended 5params: one for where photoinhibition happens
Eq12_Blackman_extended_5params <- function(data, parameters){

    # Pmax, alpha, beta, Ikb are positive parameters in the model
    Pmax   <- abs(parameters["Pmax"])
    alpha  <- abs(parameters["alpha"])
    beta   <- abs(parameters["beta"])
    Ibreak <- abs(parameters["Ibreak"])
    R  <- parameters["R"]
    Ik <- Pmax / alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <-
        alpha * I + Heaviside(I - Ik) * (Pmax - alpha * I) +
        Heaviside(I - Ibreak) * (-beta * (I - Ibreak)) + R

    return(PPhat)
}

# Eq13: Vollenweider 1958 (5params)
Eq13_Vollenweider <- function(data, parameters){

    # Pmax, alpha, and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    shape <- parameters["shape"]
    R  <- parameters["R"]

    Ik  <- Pmax/alpha
    Ikb <- Pmax/beta
    n   <- shape

    # irradiance profile
    I <- data$I

    Iratio_a <- I / alpha
    Iratio_b <- I / beta

    # calculate primary production values
    PPhat <-
        Pmax * ( Iratio_a / sqrt(1 + Iratio_a^2) ) * ( 1 / sqrt( 1 + Iratio_b^2 ) )^n + R

    return(PPhat)
}

# Eq14: Steele 1962 (3 params)
Eq14_Steele <- function(data, parameters){

    # Pmax, and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax / alpha

    # irradiance profile
    I <- data$I

    Iratio_a <- I / Ik

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a * exp(1 - Iratio_a) + R

    return(PPhat)
}

# Eq15: Steele 1962 (4 params)
Eq15_Steele_modified <- function(data, parameters){

    # Pmax, and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / Ik
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a * exp(1 - Iratio_b) + R

    return(PPhat)
}

# Eq16: Peeters 1978 (4 params)
Eq16_Peeters <- function(data, parameters){

    # Pmax, and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / Ik
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a / ( Iratio_a * Iratio_b + Iratio_a + 1 ) + R

    return(PPhat)
}

# Eq17: Platt 1980 (4 params)
Eq17_Platt <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / Ik
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * (1 - exp(-Iratio_a)) * exp(-Iratio_b) + R

    return(PPhat)
}

# Eq18: Neale 1986 (5 params)
Eq18_Neale <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax   <- abs(parameters["Pmax"])
    alpha  <- abs(parameters["alpha"])
    beta   <- abs(parameters["beta"])
    Ibreak <- abs(parameters["Ibreak"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / Ik
    Iratio_b <- (I - Ibreak) / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * exp(-Iratio_b) + R

    return(PPhat)
}

# Eq19: Neale 1986 (4 params)
Eq19_Neale_modified <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    I_ratio_1 <- I / Ik
    I_ratio_2 <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * tanh(I_ratio_1) * exp(-I_ratio_2) + R

    return(PPhat)
}

# Eq20: 2013 year (4 params)
Eq20_Ye <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- 1 + (I / Ik)
    Iratio_b <- 1 - (I / Ikb)

    # calculate primary production values
    PPhat <-
        alpha * I * ( Iratio_a / Iratio_b ) + R

    return(PPhat)
}

# Eq18: ??? year (4 params)
Eq21_Baly_extended <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / (I + Ik)
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a * exp(-Iratio_b) + R

    return(PPhat)
}

# Eq19: ??? year (4 params)
Eq22_Smith_extended <- function(data, parameters){

    # Pmax, and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / sqrt(I^2 + Ik^2)
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a * exp(-Iratio_b) + R

    return(PPhat)
}

# Eq23: ??? year (5 params)
Eq23_Prioul_extended <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    b   <- parameters["shape"]
    R   <- parameters["R"]

    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # shape param has to be between 0 to 1 based on the author
    theta <- sigmoid_fn(b)

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    Iratio_a <- (I/Ik) + 1
    Iratio_b <- I / Ikb

    # lower roots of a second order polynomial function
    I_ratio_Lroots <-  Iratio_a - sqrt(Iratio_a^2 - 4 * theta * ( Iratio_a - 1 ))

    # calculate primary production values
    PPhat <-
        Pmax * I_ratio_Lroots * exp(-Iratio_b) + R

    return(PPhat)
}

# Eq24: ??? year (5 params)
Eq24_Bannister_extended <- function(data, parameters){

    # Pmax, and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    b  <- abs(parameters["shape"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a <- I / ( (I^b + Ik^b)^(1/b) )
    Iratio_b <- I / Ikb

    # calculate primary production values
    PPhat <-
        Pmax * Iratio_a * exp(-Iratio_b) + R

    return(PPhat)
}

## 2023 4params MODELS #########------------
Eq25_Tanh_Tanh_1 <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I

    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * tanh(Iratio_b) + R

    return(PPhat)
}

Eq26_Tanh_Tanh_2 <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I

    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * tanh(Iratio_b)^2 + R

    return(PPhat)
}

# 5 params
Eq27_Tanh_Tanh_shape <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    shape <- abs(parameters["shape"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I
    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * tanh(Iratio_b)^shape + R

    return(PPhat)
}

Eq28_Tanh_Expo <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b  <- Ikb / I
    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * (1 - exp(-Iratio_b)) + R

    return(PPhat)
}

Eq29_Expo_Tanh_1 <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I
    # calculate primary production values
    PPhat <-
        Pmax *  ( 1 - exp(-Iratio_a) ) * tanh(Iratio_b) + R

    return(PPhat)
}

Eq30_Expo_Tanh_2 <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I
    # calculate primary production values
    PPhat <-
        Pmax *  ( 1 - exp(-Iratio_a) ) * tanh(Iratio_b)^2 + R

    return(PPhat)
}

# 5 params
Eq31_Expo_Tanh_shape <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    shape <- abs(parameters["shape"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I

    # calculate primary production values
    PPhat <-
        Pmax *  ( 1 - exp(-Iratio_a) ) * tanh(Iratio_b)^shape + R

    return(PPhat)
}

Eq32_Expo_Expo <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I
    # calculate primary production values
    PPhat <-
        Pmax *  ( 1 - exp(-Iratio_a) ) * ( 1 - exp(-Iratio_b) ) + R

    return(PPhat)
}

Eq33_Amirian_Tanh <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b  <- Ikb / I
    # calculate primary production values
    PPhat <-
        # Pmax * tanh(Iratio_a) * tanh(Iratio_b)^(cosh(1)) + R
        Pmax * tanh(Iratio_a) * tanh(Iratio_b)^(cosh(1)^2) + R

    return(PPhat)
}

Eq34_Amirian_Expo <- function(data, parameters){

    # Pmax, alpha and beta are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    beta  <- abs(parameters["beta"])
    R   <- parameters["R"]
    Ik  <- Pmax / alpha
    Ikb <- Pmax / beta

    # irradiance profile
    I <- data$I

    Iratio_a  <- I / Ik
    Iratio_b <- Ikb / I
    # calculate primary production values
    PPhat <-
        # Pmax * ( 1 - exp(-Iratio_a) ) * tanh(Iratio_b)^cosh(2) + R
        Pmax * ( 1 - exp(-Iratio_a) ) * tanh(Iratio_b)^(cosh(1)^2) + R

    return(PPhat)
}
