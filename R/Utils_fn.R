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
Model_setup <- function(which_model, data_type){

    # Table 1. light-limited case ----
    if(nameFormat(which_model) %in% c("linear", "lm") &
       nameFormat(data_type) %in% c("lightlimited", "lm", "ll", "llimited")){
        return(lm_piCurve)}

    # Table 2. light-saturating case ----
    if(nameFormat(data_type) %in% c("light-saturating", "saturating", "ls")){

        # LS1
        if(nameFormat(which_model) %in% c("blackmanpw", "blackman", "piecewise", "pw", "ls1", 1, "1"))
        {return(LS1_Blackman)}

        # LS2
        if(nameFormat(which_model) %in% c("balyrh", "baly", "ls2", 2, "2"))
        {return(LS2_Baly)}

        # LS3
        if(nameFormat(which_model) %in% c("smithrh", "smith", "ls3", 3, "3"))
        {return(LS3_Smith)}

        # LS4
        if(nameFormat(which_model) %in% c("webbexp", "webb", "ls4", 4, "4"))
        {return(LS4_Webb)}

        # LS5
        if(nameFormat(which_model) %in% c("jassbytanh", "jassby", "tanh", "ls5", 5, "5"))
        {return(LS5_Jassby)}

        # LS6
        if(nameFormat(which_model) %in% c("prioulnonrh", "prioul", "ls6", 6, "6"))
        {return(LS6_Prioul)}

        # LS7
        if(nameFormat(which_model) %in% c("bannistergrh", "bannister", "ls7", 7, "7"))
        {return(LS7_Bannister)}
    }

    # Table 3. photoinhibition models ----
    if(nameFormat(data_type) %in% c("photoinhibition", "ph")){
    # Ph01
    if(nameFormat(which_model) %in% c("steeleexp", "steele", "ph01", 1, "01", "1"))
    {return(Ph01_Steele)}

    # Ph02
    if(nameFormat(which_model) %in% c("peetersrational", "peeters", "ph02", 2, "02", "2"))
    {return(Ph02_Peeters)}

    # Ph03
    if(nameFormat(which_model) %in% c("plattexp", "platt", "ph03", 3, "03", "2"))
    {return(Ph03_Platt)}

    # Ph04
    if(nameFormat(which_model) %in% c("nealeexp", "neale", "ph04", 4, "04", "4"))
    {return(Ph04_Neale_modified)}

    # ph05
    if(nameFormat(which_model) %in% c("balyexp", "baly", "ph05", 5, "05", "5"))
    {return(Ph05_Baly_extended)}

    # ph06
    if(nameFormat(which_model) %in% c("smithexp", "smith", "ph06", 6, "06", "6"))
    {return(Ph06_Smith_extended)}

    # Ph07
    if(nameFormat(which_model) %in% c("bannisterexp", "bannister", "ph07", 7, "07", "7"))
    {return(Ph07_Bannister_extended)}

    # Ph08
    if(nameFormat(which_model) %in% c("prioulexp", "prioul", "ph08", 8, "08", "8"))
    {return(Ph08_Prioul_extended)}

    # Ph09
    if(nameFormat(which_model) %in% c("extendedblackman", "blackman", "pice-wise", "pw", "ph09", 9, "09", "9"))
    {return(Ph09_Blackman_extended_4params)}

    # Ph10
    if(nameFormat(which_model) %in% c("doubletanh", "2tanh", "ph10", 10, "10"))
    {return(Ph10_double_tanh)}

    # Ph11
    if(nameFormat(which_model) %in% c("doubletanhshp", "doubletanhshape",
                                      "2tanhshp", "2tanhshape", "ph11", 11, "11"))
    {return(Ph11_double_tanh_shp)}

    # Ph12
    if(nameFormat(which_model) %in% c("exptanh", "ph12", 12, "12"))
    {return(Ph12_exp_tanh)}

    # Ph13
    if(nameFormat(which_model) %in% c("exptanhshp", "exptanhshape", "ph13", 13, "13"))
    {return(Ph13_exp_tanh_shp)}

    # Ph14
    if(nameFormat(which_model) %in% c("tanhexp", "tanhrcpexp", "ph14", 14, "14"))
    {return(Ph14_tanh_rcp_exp)}

    # Ph15
    if(nameFormat(which_model) %in% c("doubleexp", "exprcpexp",  "2exp", "ph15", 15, "15"))
    {return(Ph15_exp_rcp_exp)}
        }
}

# Log likelihood function
LogL <- function(data, model_fit, StDev){
    # normal log-likelihood function
    sum(dnorm(data$PP, mean = model_fit, sd = StDev, log = TRUE))
}

add_shape_par <- function(parameters){c(parameters, shape = 1)}

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