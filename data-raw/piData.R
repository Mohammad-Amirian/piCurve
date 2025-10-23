## code to prepare `piData` dataset goes here


piDataSet <- read.csv("data-raw/piDataSet.csv")

usethis::use_data(piDataSet, overwrite = TRUE)
