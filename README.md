## Overview

The piCurve package offers a comprehensive suite of
photosynthesis–irradiance (PI) models in a user-friendly environment,
enabling users to explore how varying irradiance levels influence
photosynthetic or growth rates using statistical methods. The models
have been rigorously tested and validated against experimental data,
ensuring both reliability and accuracy. In total, 24 PI models are
included this library (`?Model_piCurve`). Optimal parameters for each
model formulation are estimated through non-linear optimization using
two statistical approaches: mean squared error (MSE) and maximum
likelihood estimation (MLE). A built-in dataset (`?piDataSet`)
containing eight independent PI incubation samples is included for model
testing and validation.

## Installation

Install the package from GitHub:

``` r
remotes::install_github("Mohammad-Amirian/piCurve")
```

    ## rlang      (1.1.6    -> 1.1.7   ) [CRAN]
    ## lifecycle  (1.0.4    -> 1.0.5   ) [CRAN]
    ## cpp11      (0.5.1    -> 0.5.2   ) [CRAN]
    ## vctrs      (0.6.5    -> 0.7.0   ) [CRAN]
    ## pillar     (1.11.0   -> 1.11.1  ) [CRAN]
    ## magrittr   (2.0.3    -> 2.0.4   ) [CRAN]
    ## timeDate   (4041.110 -> 4051.111) [CRAN]
    ## S7         (0.2.0    -> 0.2.1   ) [CRAN]
    ## isoband    (0.2.7    -> 0.3.0   ) [CRAN]
    ## tibble     (3.3.0    -> 3.3.1   ) [CRAN]
    ## gss        (2.2-9    -> 2.2-10  ) [CRAN]
    ## timeSeries (4041.111 -> 4052.112) [CRAN]
    ## ggplot2    (3.5.2    -> 4.0.1   ) [CRAN]
    ## fBasics    (4041.97  -> 4052.98 ) [CRAN]

    ## 
    ## The downloaded binary packages are in
    ##  /var/folders/fl/5wfmmnsn5d50x8xlj07dmms80000gn/T//RtmpMvpaJm/downloaded_packages
    ## ── R CMD build ─────────────────────────────────────────────────────────────────
    ##      checking for file ‘/private/var/folders/fl/5wfmmnsn5d50x8xlj07dmms80000gn/T/RtmpMvpaJm/remotes2a954d61b93f/Mohammad-Amirian-piCurve-0333299/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/fl/5wfmmnsn5d50x8xlj07dmms80000gn/T/RtmpMvpaJm/remotes2a954d61b93f/Mohammad-Amirian-piCurve-0333299/DESCRIPTION’
    ##   ─  preparing ‘piCurve’:
    ##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
    ##   ─  checking for LF line-endings in source and make files and shell scripts
    ##   ─  checking for empty or unneeded directories
    ##   ─  building ‘piCurve_0.3.3.tar.gz’
    ##      
    ## 

## How to use piCurve

``` r
library(piCurve)
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.3

### Classfying Data

To classify PI data as light-limited (ll), light-saturated (ls), or
photoinhibited (ph), you can use the `DataType_piCurve()` function.
Below, we demonstrate how to apply this function to the built-in
dataset.

``` r
# Split the data
grouped <- 
    piDataSet |> 
    dplyr::group_by(pi_number) |> 
    dplyr::group_split()

# Extract pi_number values
pi_numbers <- 
    piDataSet |> 
    dplyr::group_by(pi_number) |> 
    dplyr::group_keys() |> 
    dplyr::pull(pi_number)

# Apply function and combine with pi_number
result_df <- data.frame(
    pi_number = pi_numbers,
    data_type = unlist(parallel::mclapply(grouped, function(x) DataType_piCurve(data = x, n_cores = 1)$data_type, mc.cores = 6))
)
 
dplyr::count(result_df, data_type)
```

    ##   data_type n
    ## 1        ls 4
    ## 2        ph 4

# Model Fitting

To fit a PI model to a given data sample, you can use the
`Fit_piModel()` function. In the example below, we use the default
settings of this function, which fits the Amirian–tanh model (Ph10; see
`?Model_piCurve()`) to the data, assuming a dark reaction rate of zero
(R = 0). This model is chosen because it provides coherent and
consistent outputs. For more details, please refer to the associated
reference.

``` r
df <- grouped[[5]]

# fit model on a given data
fit <- Fit_piModel(data = df, STATapp = "MLE", Hessian = TRUE)

print(fit)
```

    ## $par
    ##        Pmax       alpha        beta       StDev 
    ## 1.852016996 0.061371599 0.003337527 0.063776649 
    ## 
    ## $value
    ## [1] -80.00742
    ## 
    ## $counts
    ## function gradient 
    ##      325       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## $hessian
    ##                Pmax        alpha         beta         StDev
    ## Pmax    7615.940697 1.139034e+04 8.674185e+05 -1.080748e+00
    ## alpha  11390.339148 5.140572e+05 5.835343e-02  1.240775e+02
    ## beta  867418.505124 5.835343e-02 1.647622e+08  1.170479e+06
    ## StDev     -1.080748 1.240775e+02 1.170479e+06  2.956539e+04
    ## 
    ## $info_criteria
    ##       AIC      AICc       BIC 
    ## -152.0148 -151.2876 -143.6375 
    ## 
    ## $model
    ## [1] "ph10"
    ## 
    ## $SQA
    ##        R2     R2adj 
    ## 0.9925895 0.9921925

### Tidy Fit

You can use the `Tidy_piCurve()` function to convert the output of your
model fit into a tidy tibble format.

``` r
# tidy the output of the model
Tidy_piCurve(fit)$params
```

    ## # A tibble: 4 × 5
    ##   params Est_Val  Est_Err `LB_95%` `UB_95%`
    ##   <chr>    <dbl>    <dbl>    <dbl>    <dbl>
    ## 1 Pmax   1.85    0.0315    1.79     1.91   
    ## 2 alpha  0.0614  0.00156   0.0583   0.0644 
    ## 3 beta   0.00334 0.000248  0.00285  0.00382
    ## 4 StDev  0.0638  0.0114    0.0414   0.0862

### Plotting

You can use the `Plot_piCurve()` function to visualize the fitted model
output along with the predicted interval (95% by default) over the data.

``` r
p <- # plot the fitted model 
    Plot_piCurve(fit, df, length_out = 200, add_CI = TRUE, n_cores = 6)

p + # increase the size of axis title and font
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

### Other Useful Functions

In addition to the core functions described above, the piCurve package
includes several other helpful tools for model analysis and
interpretation. A few of them are highlighted below.

If the user forgot to set `Hessian = TRUE` when using the
`Fit_piModel()` function, the information matrix can still be calculated
afterward using the `InfoMat_piCurve()` function, as shown below:

``` r
InfoMat_piCurve(parameters = fit$par, model_name = fit$model, data = df)
```

    ##                Pmax         alpha          beta         StDev
    ## Pmax   7.615941e+03  1.139034e+04 -8.674185e+05 -1.080748e+00
    ## alpha  1.139034e+04  5.140572e+05 -5.835343e-02  1.240775e+02
    ## beta  -8.674185e+05 -5.835343e-02  1.647622e+08 -1.170479e+06
    ## StDev -1.080748e+00  1.240775e+02 -1.170479e+06  2.956539e+04

To add confidence intervals to the predicted response variable (`y`),
you can use the `addCI_to_piPred()` function.

``` r
addCI_to_piPred(Fitted_Model = fit, irrad = df$I) |> head(3)
```

    ## # A tibble: 3 × 3
    ##   irradiance    ymin    ymax
    ##        <dbl>   <dbl>   <dbl>
    ## 1     0.0846 0.00495 0.00544
    ## 2     0.127  0.00742 0.00816
    ## 3     0.169  0.00990 0.0109

To increase the resolution of prediction value across a finer irradiance
space, you can use `?highRes_piPred()` function. List of models are also
listed in `?Model_piCurve()`.

## Phytoplankton Absorption Coefficient

Chlorophyll-a–specific phytoplankton absorption,
*a*<sub>*p**h*</sub><sup>\*</sup>(*λ*), and the corresponding
phytoplankton absorption coefficient, *a*<sub>*p**h*</sub>(*λ*), can be
estimated from chl-a concentration using the empirical spectral
parameterization of Bricaud et al. (1995). This approach provides a
globally averaged relationship between chlorophyll concentration and
phytoplankton absorption in the absence of direct optical measurements
(Details `?absCoef`). Below is an example illustrating the computation
and visualization of *a*<sub>*p**h*</sub>(*λ*) for a chl-a concentration
of 0.9 mg m<sup>−3</sup>.

``` r
absCoef(chla = 0.9) |>
    ggplot(aes(x = lambda_nm, y = a_ph)) +
    geom_line() +
    theme_bw(base_size = 14) +
    labs(
        x = "Wavelength (nm)",
        y = expression(a[ph](lambda)~~(m^{-1}))
    )
```

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

## Citation

This package is developed based on the research outlined in the
following paper. If you utilize the `piCurve` package in your work, we
kindly request that you cite the referenced publication.

M. Amirian, M., Finkel, Z.V., Devred, E., Irwin, A.J. “*Parameterization
of photoinhibition for phytoplankton*”. Commun Earth Environ 6, 707
(2025). <https://doi.org/10.1038/s43247-025-02686-3>

M. Amirian, M. and Irwin, A.J. “*piCurve: an R package for modeling
photosynthesis–irradiance curves*”. arXiv preprint arXiv:2508.14321v1,
(2025). <https://arxiv.org/abs/2508.14321>
