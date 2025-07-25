## Overview

The piCurve package offers a comprehensive suite of
photosynthesis–irradiance (PI) models in a user-friendly environment,
enabling users to explore how varying irradiance levels influence
photosynthetic or growth rates using statistical methods. The models
have been rigorously tested and validated against experimental data,
ensuring both reliability and accuracy. Optimal parameters for each
model formulation are estimated through non-linear optimization using
two statistical approaches: mean squared error (MSE) and maximum
likelihood estimation (MLE). A built-in dataset (`?piDataSet`)
containing eight independent PI incubation samples is included for model
testing and validation.

## How to use piCurve

``` r
rm(list = ls())
devtools::load_all()
```

    ## ℹ Loading piCurve

``` r
library(piCurve)
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.3.3

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(fBasics)
```

    ## Warning: package 'fBasics' was built under R version 4.3.3

### Classfying Data

To classify PI data as light-limited (ll), light-saturated (ls), or
photoinhibited (ph), you can use the `DataType_piCurve()` function.
Below, we demonstrate how to apply this function to the built-in
dataset.

``` r
# Split the data
grouped <- 
    piDataSet |> 
    group_by(pi_number) |> 
    group_split()

# Extract pi_number values
pi_numbers <- 
    piDataSet |> 
    group_by(pi_number) |> 
    group_keys() |> 
    pull(pi_number)

# Apply function and combine with pi_number
result_df <- data.frame(
    pi_number = pi_numbers,
    data_type = sapply(grouped, function(x) DataType_piCurve(data = x, n_cores = 6)$data_type)
)
 
count(result_df, data_type)
```

    ##   data_type n
    ## 1        ls 5
    ## 2        ph 3

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
    ## 1.851955156 0.061374091 0.003337316 0.063783954 
    ## 
    ## $value
    ## [1] -80.00741
    ## 
    ## $counts
    ## function gradient 
    ##      396       NA 
    ## 
    ## $convergence
    ## [1] 0
    ## 
    ## $message
    ## NULL
    ## 
    ## $hessian
    ##               Pmax        alpha         beta        StDev
    ## Pmax  7.614233e+03 1.138670e+04 8.671789e+05 4.957483e+00
    ## alpha 1.138670e+04 5.138762e+05 5.825017e-02 1.059650e+02
    ## beta  8.671789e+05 5.825017e-02 1.647233e+08 1.170477e+06
    ## StDev 4.957483e+00 1.059650e+02 1.170477e+06 2.954846e+04
    ## 
    ## $info_criteria
    ##       AIC      AICc       BIC 
    ## -150.0148 -148.9037 -139.5431 
    ## 
    ## $model
    ## [1] "2tanh"
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

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

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
    ## Pmax   7.614233e+03  1.138670e+04 -8.671789e+05  4.957483e+00
    ## alpha  1.138670e+04  5.138762e+05 -5.825017e-02  1.059650e+02
    ## beta  -8.671789e+05 -5.825017e-02  1.647233e+08 -1.170477e+06
    ## StDev  4.957483e+00  1.059650e+02 -1.170477e+06  2.954846e+04

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

## Citation

This package is developed based on the research outlined in the
following paper. If you utilize the `piCurve` package in your work, we
kindly request that you cite the referenced publication.

Mohammad Amirian, Emmanuel Devred, Zoe V. Finkel, Andrew J. Irwin. “*A
new parameterization of photoinhibition for phytoplankton*,” arXiv
(2024) 1–33. 10.48550/arXiv.2412.17923.
