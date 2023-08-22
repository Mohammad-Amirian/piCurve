#####################################################################
# Functions used to model pi-curve with no photoinhibition (beta = 0)

################## model notations ##################################
# Pmax : maximum photosynthesis rate (sign +)
# alpha: light-saturation slope, low light level (sign +)
# R    : dark reaction parameter (sign + OR -)
# b    : shape parameter (sign +)
#####################################################################

# Eq1: linear regression
Eq1_lm <- function(data, parameters){

    # alpha are positive parameters in the model
    alpha <- parameters["alpha"]
    R  <- parameters["R"]

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <- alpha * I + R

    return(PPhat)
}

# Eq1: Blackman 1905 (Bi-linear)
Eq2_Blackman <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <- alpha*I + Heaviside(I - Ik) * (Pmax - alpha*I) + R

    return(PPhat)
    }

# Eq2: Baly 1935 (rectangular hyperbola)
Eq3_Baly <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / (I + Ik)

    # calculate primary production values
    PPhat <- Pmax * I_ratio + R

    return(PPhat)
}

# Eq3: Smith 1936 (modified rectangular hyperbola)
Eq4_Smith <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / sqrt(I^2 + Ik^2)

    # calculate primary production values
    PPhat <- Pmax * I_ratio + R

    return(PPhat)
}

# Eq5: Talling 1957 (modified--original suggestion belongs to Vollenweider 1958)
Eq5_Talling <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / Ik

    # calculate primary production values
    PPhat <- Pmax * log(2 * I_ratio + 1) + R

    return(PPhat)
}

# Eq5: Vollenweider 1958
Eq6_Vollenweider_Ln <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / Ik

    # calculate primary production values
    PPhat <- Pmax * log( I_ratio + sqrt(1 + I_ratio^2) ) + R

    return(PPhat)
}

# Eq4: Webb 1974 (exponential)
Eq7_Webb <- function(data, parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <- Pmax * ( 1 - exp (- I/Ik ) ) + R

    return(PPhat)
}

# Eq5: Jassby 1976 (hyperbolic tangent)
Eq8_Jassby <- function(data,  parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate primary production values
    PPhat <- Pmax * tanh( I/Ik ) + R

    return(PPhat)
}

# Eq6: Prioul 1977 (non-rectangular hyperbola)
Eq9_Prioul <- function(data,  parameters){

    # Pmax and alpha are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    R  <- parameters["R"]
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
    PPhat <- Pmax * I_ratio_Lroots / (2 * theta) + R

    return(PPhat)
}

# Eq7: Bannister 1979 (generalized rectangular hyperbola)
Eq10_Bannister <- function(data, parameters){

    # Pmax, alpha, and shape are positive parameters in the model
    Pmax  <- abs(parameters["Pmax"])
    alpha <- abs(parameters["alpha"])
    b  <- abs(parameters["shape"])
    R  <- parameters["R"]
    Ik <- Pmax/alpha

    # irradiance profile
    I <- data$I

    # calculate irradiance ratio
    I_ratio <- I / ( (I^b + Ik^b)^(1/b) )

    # calculate primary production values
    PPhat <- Pmax * I_ratio + R

    return(PPhat)
}
