.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        # "piCurve 0.2.4 loaded:: Windows parallel bug in Fit_piModel() fixed."
        "piCurve 0.2.5 loaded:: DataType_piCurve handles imbalanced more delicately now using both adjR2 & AIC"
    )
}
