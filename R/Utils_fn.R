# Utility functions

# Sigmoid function, mapping input values to a range between 0 and 1
sigmoid_fn <- function(x) {
    return(1 / (1 + exp(-x)))
}

# remove spacing and dash/under line(s) from the name
nameFormat <- function(model_name){
    x <- gsub(" ", "", model_name) # rm spacing
    x <- gsub("_", "", x)          # rm underline
    x <- gsub("-", "", x)          # rm dash line
    return(tolower(x))
}

# which equation?
Model_setup <- function(which_model){


    # No photo-inhibition case
    # Eq1 ----
    if(nameFormat(which_model) %in% c("eq1lm", "lmeq1", "eq1", "lm", 1, "1"))
    {return(Eq1_lm)}

    # Eq2 ----
    if(nameFormat(which_model) %in% c("eq2blackman", "blackmaneq2", "eq2", 2, "2"))
    {return(Eq2_Blackman)}

    # Eq3 ----
    if(nameFormat(which_model) %in% c("eq3baly", "balyeq3", "eq3", 3, "3"))
    {return(Eq3_Baly)}

    # Eq4 ----
    if(nameFormat(which_model) %in% c("eq4smith", "smitheq4", "eq4", 4, "4"))
    {return(Eq4_Smith)}

    # Eq5  ----
    if(nameFormat(which_model) %in% c("eq5talling", "tallingeq5", "eq5", 5, "5"))
    {return(Eq5_Talling)}

    # Eq6 ----
    if(nameFormat(which_model) %in% c("eq6vollenweiderln", "eq6vollenweider",
                                   "vollenweiderlneq6", "vollenweidereq6",
                                   "Eq6", 6, "6"))
    {return(Eq6_Vollenweider_Ln)}

    # Eq7 ----
    if(nameFormat(which_model) %in% c("eq7webb", "webbeq7", "eq7", 7, "7"))
    {return(Eq7_Webb)}

    # Eq8 ----
    if(nameFormat(which_model) %in% c("eq8jassby", "jassbyeq8", "eq8", 8, "8"))
    {return(Eq8_Jassby)}

    # Eq9 ----
    if(nameFormat(which_model) %in% c("eq9prioul", "priouleq9", "eq9", 9, "9"))
    {return(Eq9_Prioul)}

    # Eq10 ----
    if(nameFormat(which_model) %in% c("eq10bannister", "bannistereq10", "eq10", 10, "10"))
    {return(Eq10_Bannister)}


    # photoinihibtion models ----
    # Eq11 ----
    if(nameFormat(which_model) %in% c("eq11blackman", "blackmaneq11",
                                      "eq11blackmanextended1", "blackmanextended1eq11",
                                      "blackmanextended1", "extended1blackman",
                                      "eq11", 11, "11"))
    {return(Eq11_Blackman_extended_4params)}

    # Eq12 ----
    if(nameFormat(which_model) %in% c("eq12blackman", "blackmaneq12",
                                   "eq12blackmanextended2", "blackmanextended2eq12",
                                   "blackmanextended2", "extended2blackman",
                                   "eq12", 12, "12"))
    {return(Eq12_Blackman_extended_5params)}

    # Eq13 ----
    if(nameFormat(which_model) %in% c("eq13vollenweider", "vollenweidereq13",
                                   "eq13", 13, "13"))
    {return(Eq13_Vollenweider)}

    # Eq14 ----
    if(nameFormat(which_model) %in% c("eq14steele", "steeleeq14",
                                   "eq14", 14, "14"))
    {return(Eq14_Steele)}

    # Eq15 ----
    if(nameFormat(which_model) %in% c("eq15steelemodified", "steelemodifiedeq15",
                                      "eq15", 15, "15"))
    {return(Eq15_Steele_modified)}

    # Eq16 ----
    if(nameFormat(which_model) %in% c("eq16peeters", "peeterseq16",
                                      "eq16", 16, "16"))
    {return(Eq16_Peeters)}

    # Eq17 ----
    if(nameFormat(which_model) %in% c("eq17platt", "platteq17",
                                      "eq17", 17, "17"))
    {return(Eq17_Platt)}

    # Eq18 ----
    if(nameFormat(which_model) %in% c("eq18neale", "nealeeq18",
                                      "eq18", 18, "18"))
    {return(Eq18_Neale)}
    # Eq19 ----
    if(nameFormat(which_model) %in% c("eq19nealemodified", "nealemodifiedeq19",
                                      "eq19", 19, "19"))
    {return(Eq19_Neale_modified)}
    # Eq20 ----
    if(nameFormat(which_model) %in% c("eq20ye", "yeeq20",
                                      "eq20", 20, "20"))
    {return(Eq20_Ye)}
    # Eq21 ----
    if(nameFormat(which_model) %in% c("eq21balyextended", "balyextendedeq21",
                                      "eq21", 21, "21"))
    {return(Eq21_Baly_extended)}

    # Eq22 ----
    if(nameFormat(which_model) %in% c("eq22smithextended", "smithextendedeq22",
                                      "eq22", 22, "22"))
    {return(Eq22_Smith_extended)}

    # Eq23 ----
    if(nameFormat(which_model) %in% c("eq23prioulextended", "prioulextendedeq23",
                                      "eq23", 23, "23"))
    {return(Eq23_Prioul_extended)}

    # Eq24 ----
    if(nameFormat(which_model) %in% c("eq24bannisterextended", "bannisterextendedeq24",
                                      "eq24", 24, "24"))
    {return(Eq24_Bannister_extended)}

    # Eq25 ----
    if(nameFormat(which_model) %in% c("eq25tanhtanh1", "tanhtanh1eq25",
                                      "eq25", 25, "25"))
    {return(Eq25_Tanh_Tanh_1)}

    # Eq26 ----
    if(nameFormat(which_model) %in% c("eq26tanhtanh2", "tanhtanh2eq26",
                                      "eq26", 26, "26"))
    {return(Eq26_Tanh_Tanh_2)}

    # Eq27 ----
    if(nameFormat(which_model) %in% c("eq27tanhtanhshape", "tanhtanhshapeeq27",
                                      "eq27", 27, "27"))
    {return(Eq27_Tanh_Tanh_shape)}

    # Eq28 ----
    if(nameFormat(which_model) %in% c("eq28tanhexpo", "tanhexpoeq28",
                                      "eq28", 28, "28"))
    {return(Eq28_Tanh_Expo)}

    # Eq29 ----
    if(nameFormat(which_model) %in% c("eq29expotanh1", "expotanh1eq29",
                                      "eq29", 29, "29"))
    {return(Eq29_Expo_Tanh_1)}

    # Eq30 ----
    if(nameFormat(which_model) %in% c("eq30expotanh2", "expotanh2eq30",
                                      "eq30", 30, "30"))
    {return(Eq30_Expo_Tanh_2)}

    # Eq31 ----
    if(nameFormat(which_model) %in% c("eq31expotanhshape", "expotanhshapeeq31",
                                      "eq31", 31, "31"))
    {return(Eq31_Expo_Tanh_shape)}

    # Eq32 ----
    if(nameFormat(which_model) %in% c("eq32expoexpo", "expoexpoeq32",
                                      "eq32", 32, "32"))
    {return(Eq32_Expo_Expo)}

    # Eq33 ----
    if(nameFormat(which_model) %in% c("eq33amiriantanh", "amiriantanheq33",
                                      "eq33", 33, "33"))
    {return(Eq33_Amirian_Tanh)}

    # Eq34 ----
    if(nameFormat(which_model) %in% c("eq34amirianexpo", "amirianexpoeq34",
                                      "eq34", 34, "34"))
    {return(Eq34_Amirian_Expo)}
}

# Log likelihood function
LogL <- function(data, model_fit, StDev){
    # normal log-likelihood function
    sum(dnorm(data$PP, mean = model_fit, sd = StDev, log = TRUE))
}

add_shape_par <- function(parameters){c(parameters, shape = 1)}

# Eq12_Blackman_extended_5params and Eq18_Neale have an Ibreak parameter, where
# photoinhibition kicks. Wright chose of initial input in optimization is crucial.
# This function find a proper starting point using Eq11
Approx_Ibreak_par <- function(parameters, data){
    fit <- OptPIparams(parameters, model_name = "eq11", data)
    Ibreak <- as.numeric(fit$par["Pmax"] / fit$par["beta"])

    return(c(fit$par, Ibreak = Ibreak))
}

# reshape the fit output
reshape_mse <-
    function(fitted_model){
        c(Model = fitted_model$model,fitted_model$par, MSE = fitted_model$value,
          fitted_model$SQA, Convrg = fitted_model$convergence)}

reshape_piFits <-
    function(fitted_model, STATapp) {
        ifelse(
            missing(STATapp) || tolower(STATapp) == "mse",
            out <- reshape_mse(fitted_model),
            out <- c(reshape_mse(fitted_model), fitted_model$info_criteria)
        )
        return(out)
    }
