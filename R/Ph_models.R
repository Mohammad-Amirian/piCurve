#####################################################################
# Functions used to model pi-curve with with photoinhibition (beta ~= 0)
# List of models are given in '?piCurve::Model_piCurve'

################## model notations ##################################
# Pmax : maximum photosynthesis rate (sign +)
# alpha: slope of the light-saturation curve at low light (sign +)
# beta : photoinhibition parameter (sign +)
# R    : dark reaction rate (sign + OR -)
# shape: shape parameter (sign +)

#####################################################################
################## formulations #####################################

# Ph01: Steele 1962 (3 params)
Ph01_Steele <- function(data, parameters){

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

# Ph02: Peeters 1978 (4 params)
Ph02_Peeters <- function(data, parameters){

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

# Ph03: Platt 1980 (4 params)
Ph03_Platt <- function(data, parameters){

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

# Ph04: Neale 1986 (4 params)
Ph04_Neale_modified <- function(data, parameters){

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

# Ph05: ??? year (4 params)
Ph05_Baly_extended <- function(data, parameters){

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

# Ph06: ??? year (4 params)
Ph06_Smith_extended <- function(data, parameters){

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

# Ph07: ??? year (5 params)
Ph07_Bannister_extended <- function(data, parameters){

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

# Ph08: ??? year (5 params)
Ph08_Prioul_extended <- function(data, parameters){

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

# Ph09: Blackman model, extended 4 params
Ph09_Blackman_extended_4params <- function(data, parameters){

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

# Ph10: Amirian 2024 (4 params)
Ph10_double_tanh <- function(data, parameters){

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
    Iratio_bb  <- (Iratio_b)^(cosh(1)^2)
    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * tanh(Iratio_bb) + R

    return(PPhat)
}

# Ph11: Amirian 2024 (5 params)
Ph11_double_tanh_shp <- function(data, parameters){

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
    Iratio_bb <- Iratio_b^shape

    # calculate primary production values
    PPhat <-
        Pmax * tanh(Iratio_a) * tanh(Iratio_bb) + R

    return(PPhat)
}

# Ph12: Amirian 2024 (4 params)
Ph12_exp_tanh <- function(data, parameters){

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
    Iratio_bb <- Iratio_b^(cosh(1)^2)

    # calculate primary production values
    PPhat <-
        Pmax * ( 1 - exp(-Iratio_a) ) * tanh(Iratio_bb) + R

    return(PPhat)
}

# Ph13: Amirian 2024 (5 params)
Ph13_exp_tanh_shp <- function(data, parameters){

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
    Iratio_bb  <- Iratio_b^shape

    # calculate primary production values
    PPhat <-
        Pmax *  ( 1 - exp(-Iratio_a) ) * tanh(Iratio_bb) + R

    return(PPhat)
}

# Ph14: Amirian 2024 (4 params)
Ph14_tanh_rcp_exp <- function(data, parameters){

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

# Ph15: Amirian 2024 (4 params)
Ph15_exp_rcp_exp <- function(data, parameters){

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

