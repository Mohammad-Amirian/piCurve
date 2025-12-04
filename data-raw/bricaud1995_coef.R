## data-raw/bricaud1995_coef.R
## One-time script to build internal package data

library(readr)
library(dplyr)
library(usethis)

# Path to raw CSV
infile <- "data-raw/Bricaud_1995_data.csv"

# Read and sanitize
bricaud1995_coef <- read_csv(infile, show_col_types = FALSE) |>
    rename(
        lambda_nm = lambda_nm,  # wavelength (nm)
        A = A,                  # Bricaud A(lambda)
        B = B                   # Bricaud B(lambda)
    ) |>
    mutate(
        lambda_nm = as.numeric(lambda_nm),
        A = as.numeric(A),
        B = as.numeric(B)
    ) |>
    arrange(lambda_nm)

# Basic integrity checks
stopifnot(
    all(c("lambda_nm", "A", "B") %in% names(bricaud1995_coef)),
    all(is.finite(bricaud1995_coef$lambda_nm)),
    all(bricaud1995_coef$lambda_nm > 0),
    !anyNA(bricaud1995_coef)
)

# Save as internal package dataset (.rda in data/)
usethis::use_data(bricaud1995_coef, overwrite = TRUE)
