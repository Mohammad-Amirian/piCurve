

GreggCarder1990 <- read.csv(
    "data-raw/GreggCarder1990.csv",
    stringsAsFactors = FALSE
)

# optional: convert to tibble if you want
GreggCarder1990 <- as.data.frame(GreggCarder1990)

# basic checks
stopifnot(
    all(c("lambda", "H_o", "a_oz", "a_w", "a_o") %in% names(GreggCarder1990)),
    is.numeric(GreggCarder1990$lambda),
    is.numeric(GreggCarder1990$H_o),
    is.numeric(GreggCarder1990$a_oz),
    is.numeric(GreggCarder1990$a_w),
    is.numeric(GreggCarder1990$a_o)
)

# ordering
GreggCarder1990 <- GreggCarder1990[order(GreggCarder1990$lambda), ]

usethis::use_data(GreggCarder1990, overwrite = TRUE)
