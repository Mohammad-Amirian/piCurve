#####################################################################
# Functions used to model pi-curve with no photoinhibition (beta = 0).
# List of models are given in '?piCurve::Model_piCurve'

################## model notations ##################################
# Pmax : maximum photosynthesis rate (sign +)
# alpha: slope of the light-saturation curve at low light (sign +)
# R    : dark reaction rate (sign + OR -)
# b    : shape parameter (sign +)
#####################################################################
################## formulations #####################################


# lm: linear regression
lm_piCurve <- function(data, parameters){

    # alpha are positive parameters in the model
    alpha <- parameters["alpha"]


    # irradiance profile
    I <- data$I

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- alpha * I + R
    } else {
        PPhat <- alpha * I
    }

    return(PPhat)
}

# LS1: Blackman 1905 (Bi-linear)
LS1_Blackman <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])

    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- alpha*I + Heaviside(I - Ik) * (Pmax - alpha*I) + R
    } else {
        PPhat <- alpha*I + Heaviside(I - Ik) * (Pmax - alpha*I)
    }

    return(PPhat)
    }

# LS2: Baly 1935 (rectangular hyperbola)
LS2_Baly <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])

    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / (I + Ik)

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * I_ratio + R
    } else {
        PPhat <- Pmax * I_ratio
    }


    return(PPhat)
}

# LS3: Smith 1936 (modified rectangular hyperbola)
LS3_Smith <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])

    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / sqrt(I^2 + Ik^2)

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * I_ratio + R
    } else {
        PPhat <- Pmax * I_ratio
    }

    return(PPhat)
}

# LS4: Webb 1974 (exponential)
LS4_Webb <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])

    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * ( 1 - exp (- I/Ik ) ) + R
    } else {
        PPhat <- Pmax * ( 1 - exp (- I/Ik ) )
    }



    return(PPhat)
}

# LS5: Jassby 1976 (hyperbolic tangent)
LS5_Jassby <- function(data,  parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])

    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * tanh( I/Ik ) + R
    } else {
        PPhat <- Pmax * tanh( I/Ik )
    }

    return(PPhat)
}

# Eq6: Prioul 1977 (non-rectangular hyperbola)
LS6_Prioul <- function(data,  parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    b  <- parameters["shape"]
    Ik <- Pmax/alpha

    # shape param has to be between 0 to 1 based on the author
    theta <- sigmoid_fn(b)

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- (I/Ik) + 1

    # lower roots of a second order polynomial function
    I_ratio_Lroots <-  I_ratio - sqrt(I_ratio^2 - 4 * theta * ( I_ratio - 1 ))

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * I_ratio_Lroots / (2 * theta) + R
    } else {
        PPhat <- Pmax * I_ratio_Lroots / (2 * theta)
    }


    return(PPhat)
}

# Eq7: Bannister 1979 (generalized rectangular hyperbola)
LS7_Bannister <- function(data, parameters){

    # Pmax, alpha, and shape are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    b  <- abs(parameters["shape"])
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / ( (I^b + Ik^b)^(1/b) )

    # calculate primary production values
    if ("R" %in% names(parameters)){
        R  <- parameters["R"]
        PPhat <- Pmax * I_ratio + R
    } else {
        PPhat <- Pmax * I_ratio
    }

    return(PPhat)
}
