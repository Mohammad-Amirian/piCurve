## code to prepare `piData` dataset goes here


pi_data <- read.csv("data-raw/piDataSet.csv")

usethis::use_data(pi_data, overwrite = TRUE)
