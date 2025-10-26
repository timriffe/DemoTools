#' Download or Process Sex Ratio at Birth (SRB)
#'
#' Retrieves country-specific Sex Ratio at Birth (SRB) data from the most recent installed \code{wpp} package
#' (e.g., \code{wpp2022}) or returns user-supplied SRB values for a specified set of years. If the requested
#' location is not found in WPP data, a default SRB value (1.0486) is used. Missing years are filled with the
#' mean of available SRB values.
#'
#' @param SRB Optional numeric value or vector specifying SRB(s) to use. If provided, this value is repeated for all \code{DatesOut}.
#' @param location Character or numeric code specifying the country or region for which SRB data should be downloaded.
#' @param DatesOut Numeric vector of years for which SRB values are required. Must contain at least one element.
#' @param verbose Logical; if \code{TRUE}, prints progress and diagnostic messages. Default is \code{TRUE}.
#'
#' @details
#' The function identifies the latest installed \code{wpp} package (e.g., \code{wpp2022}) and extracts single-year
#' SRB data (\code{sexRatio1}). It supports both country names and ISO numeric codes as \code{location} inputs.
#' 
#' If the specified location is not found in the WPP data, a default SRB value corresponding to a sex ratio of
#' 105 male births per 100 female births (\code{1.0486}) is returned. If certain requested years are missing from
#' the WPP dataset, the missing SRB values are imputed using the mean of available SRB values.
#'
#' @return
#' A named numeric vector of SRB values for each year in \code{DatesOut}.
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom magrittr %>% 
#'
#' @examples
#' \dontrun{
#' # Example 1: Download SRB data for France
#' srb_data <- download_SRB(location = "France",
#'                          DatesOut = c(2000, 2005, 2010))
#'
#' # Example 2: Use user-supplied constant SRB value
#' srb_data <- download_SRB(SRB = 1.05,
#'                          location = "France",
#'                          DatesOut = c(2000, 2005, 2010))
#' }
#'
#' @export

download_SRB <- function(SRB = NULL, 
                         location, 
                         DatesOut, 
                         verbose = TRUE) {
  
  SRB_default <- round((1 - 0.4886) / 0.4886, 3)
  
  # Check DatesOut
  if(length(DatesOut) < 1) {
    
    stop("DatesOut must contain at least one date.")
    
  }  
  
  # If SRB provided directly
  # ------------------------------------------------------------------ #
  # CASE 1: User-supplied SRB values (no download required)
  # ------------------------------------------------------------------ #
  if (!is.null(SRB)) {
    
    SRB <- setNames(rep(SRB, length.out = length(DatesOut)), DatesOut)
    
    return(SRB)
    
  }
  
  # ------------------------------------------------------------------ #
  # CASE 2: Retrieve SRB from latest installed WPP package
  # ------------------------------------------------------------------ #
  if(is.null(SRB)) { 
  
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  if(length(installed_wpp) == 0) stop("No WPP package installed.")
  
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  
  if(latest_wpp < 2022) {
    
    stop("No single years SRB is availabe in wpp versions earlier that wpp2022.
           Please update the wpp package.")
    
  }
  
  # Load SRB data from the selected WPP package
  data("sexRatio1", package = latest_wpp)
  
  # Determine country code (numeric or character location input)
  if(is.numeric(location)) {
    
    location_code <- location
    
  } else {
    
    location_code <- sexRatio1 %>%
      as_tibble() %>%
      filter(.data$name %in% location) %>%
      select("country_code") %>%
      unique() %>%
      as.numeric()
  }
  
  # Handle missing or invalid locations gracefully
  if(is.na(location_code)) {
    
    if(verbose) {  
      
      cat(location, "not available in wpp. Using default SRB:", SRB_default, "\n")
      
    }
    
    return(setNames(rep(SRB_default, length(DatesOut)), DatesOut))
  }
  
  # Notify user of download progress
  if(verbose) {
    cat(paste0("\nDownloading SRB for ", location, 
               " for years ", round(min(DatesOut),1), " to ", round(max(DatesOut),1), "\n"))
  }
  
  # ------------------------------------------------------------------ #
  # Extract and format SRB data for the specified country and years
  # ------------------------------------------------------------------ #
  dt <- sexRatio1 %>%
    as_tibble() %>%
    filter(.data$country_code == location_code) %>%
    pivot_longer(-c("country_code", "name"), 
                 names_to  = "year", 
                 values_to = "SRB") %>%
    mutate(year = as.numeric(.data$year)) %>%
    filter(.data$year %in% floor(DatesOut))
  
  # Vectorize SRB values by year
  years_srb <- dt$year
  SRB       <- setNames(dt$SRB, years_srb)
  
  # ------------------------------------------------------------------ #
  # Handle missing years: fill with mean SRB of available years
  # ------------------------------------------------------------------ #  
  DatesOut_floor <- floor(DatesOut)
  yrs_present    <- DatesOut_floor %in% years_srb
  
  if(any(!yrs_present)) {
    mean_srb <- mean(SRB[as.character(DatesOut_floor[yrs_present])])
    SRB      <- c(SRB, setNames(rep(mean_srb, sum(!yrs_present)), DatesOut_floor[!yrs_present]))
  }
  
  # Sort SRB by year and assign final names
  SRB <- SRB[order(as.numeric(names(SRB)))]
  SRB <- setNames(SRB, DatesOut)
  }
  
  return(SRB)
}



# lets try SRB 
# SRB <- download_SRB(SRB      = NULL,
#                     location = "Argentina",
#                     DatesOut = 1950:1952, 
#                     verbose  = TRUE)

# another setup
# SRB2 <- download_SRB(SRB      = NULL,
#                      location = "Spain",
#                      DatesOut = 2000:2030, 
#                      verbose  = TRUE)

# lets try provide SRB ourselves. All good
# SRB_user <- download_SRB(SRB = SRB[1],
#                          DatesOut = 1950)


