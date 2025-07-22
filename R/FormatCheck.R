#' Format Check for Input PI Data
#'
#' This function ensures that input data are properly structured for P–I modeling.
#' Specifically, the function checks for the presence of `P` (photosynthesis rate) and `I`
#' (irradiance) columns, removes any missing values, and coerces both variables to numeric type
#' if needed.
#'
#'
#' @param data A data frame containing the input data.
#'
#' @return A cleaned and formatted data frame named `prep_data` with numeric columns `P` and `I`, and no missing values.
#' The function removes spacing and replaces commas (,), colons (:), and semicolons (;) to dot symbol.
#'
#' @examples
#' \dontrun{
#' df <- read.csv("mydata.csv")
#' clean_df <- FormatCheck_piCurve(df)
#' }
#'
#' @export
FormatCheck_piCurve <- function(data) {
    # Check if P and I columns are present
    if (!all(c("P", "I") %in% colnames(data))) {
        message("Column names 'P' and/or 'I' not found in the dataset.")
        print("Available columns are:")
        print(colnames(data))

        p_col <- readline(prompt = "Please enter the name of the column to be used as P: ")
        i_col <- readline(prompt = "Please enter the name of the column to be used as I: ")

        # Rename columns
        if (!p_col %in% colnames(data) | !i_col %in% colnames(data)) {
            stop("The specified columns were not found in the data. Please check and try again.")
        }

        colnames(data)[colnames(data) == p_col] <- "P"
        colnames(data)[colnames(data) == i_col] <- "I"
    }

    # Remove rows with NA values
    if (anyNA(data)) {
        message("Missing values found. Removing rows with NA...")
        data <- na.omit(data)
    }

    # rm spaces, commas, colons, and semicolons
    data$P <- gsub(" ", "", data$P)
    data$I <- gsub(" ", "", data$I)

    data$P <- gsub("[,;:]", ".", data$P)
    data$I <- gsub("[,;:]", ".", data$I)

    # Convert P and I to numeric (if not already)
    data$P <- suppressWarnings(as.numeric(data$P))
    data$I <- suppressWarnings(as.numeric(data$I))

    # Check for conversion issues
    if (any(is.na(data$P)) || any(is.na(data$I))) {
        warning("Some values in P or I could not be converted to numeric and were removed.")
        data <- data[!is.na(data$P) & !is.na(data$I), ]
    }

    # Return cleaned and formatted data
    prep_data <- data
    return(prep_data)
}
