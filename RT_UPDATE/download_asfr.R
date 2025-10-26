#' Download or Process Age-Specific Fertility Rates (ASFR)
#'
#' This function returns user-supplied ASFR data formatted to a specific age
#' structure or downloads and constructs ASFR values from the most recent
#' available \code{wpp} package for a given location. It supports interpolation
#' across years and aggregation to 5-year age groups.
#'
#' When \code{Asfrmat} is supplied the function only reformats and labels age
#' rows according to \code{output}. When \code{Asfrmat} is not supplied, the
#' function:
#' \enumerate{
#'   \item Finds the latest installed \code{wpp} package (must be \eqn{\ge} 2022 for single-year data),
#'   \item Loads \code{percentASFR1dt} and \code{tfr1},
#'   \item Constructs ASFR either as single-year rates or aggregated 5-year rates,
#'   \item Interpolates across requested years \code{AsfrDatesIn} using \code{interp()}.
#' }
#'
#' @param Asfrmat Optional matrix or data frame of user-supplied ASFR values.
#'   If provided, the function returns the formatted ASFR data according to
#'   \code{output}. \strong{Expected convention: rows = ages, columns = years.}
#' @param location Character or numeric code specifying the country or region
#'   for which to download data. Required if \code{Asfrmat} is not supplied.
#' @param AsfrDatesIn Numeric vector of reference years for which ASFR data are
#'   to be downloaded or aligned.
#' @param method Character specifying interpolation method to use when ASFR data
#'   are interpolated across years. Defaults to \code{"linear"}.
#' @param output Character specifying desired output age structure:
#'   \code{"single"} (single-year ages) or \code{"5-year"} (default).
#' @param ... Additional arguments passed to the interpolation routine.
#'
#' @details
#' The function expects WPP single-age percentage ASFR data (\code{percentASFR1dt})
#' and total fertility rate (\code{tfr1}) to be available in the latest installed
#' \code{wpp} package. For the \code{"single"} option, the function constructs
#' single-year ASFR by multiplying percentage ASFR by TFR and interpolating to
#' requested years. For the \code{"5-year"} option, the function first aggregates
#' percentages to 5-year groups and then multiplies by TFR; it then interpolates
#' group-level rates to requested years. When \code{Asfrmat} is provided by the
#' user, it is converted to a matrix and rownames are set according to the
#' requested \code{output}.
#'
#' @return
#' A numeric matrix of ASFR values by age (rows) and year (columns).
#'
#' @importFrom dplyr filter select mutate group_by summarise ungroup left_join
#' @importFrom dplyr group_nest
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' #importFrom magrittr %>% 
#'
#' @examples
#' \dontrun{
#' # Download 5-year ASFR for France for three years
#' asfr_data <- download_Asfr(location = "France",
#'                            AsfrDatesIn = c(2000, 2005, 2010),
#'                            output = "5-year")
#'
#' # Use user-supplied ASFR matrix (rows = ages, cols = years)
#' Asfrmat <- matrix(runif(25), nrow = 5)
#' download_Asfr(Asfrmat = Asfrmat,
#'               AsfrDatesIn = c(2000, 2005, 2010, 2015, 2020),
#'               output = "single")
#' }
#'
#' @export


download_Asfr <- function(Asfrmat     = NULL,
                          location    = NULL,
                          AsfrDatesIn = NULL,
                          method      = "linear",
                          output      = "5-year",
                          ...) {
  
  verbose <- getOption("basepop_verbose", TRUE)
  
  # -------------------------------
  # CASE 1: user-supplied Asfrmat
  # -------------------------------
  # If the user provides Asfrmat we simply coerce to matrix, assign column names
  # (years) and set rownames (ages) according to requested output format.
  #
  # NOTE: Asfrmat is expected to be in conventional reproductive ages
  # (usually 15-49 single-year ages or 5-year groups 15-19,...). If input
  # has a different age basis, it should be transformed first.
  if(!is.null(Asfrmat)) {
    
    Asfrmat           <- as.matrix(Asfrmat)
    colnames(Asfrmat) <- AsfrDatesIn
    n                 <- nrow(Asfrmat)
    
    # Assign age labels depending on output. Keep behaviour consistent with earlier code.
    if(output == "single") {
      
      # Single-year ages: assume row 1 corresponds to age 10 (user responsibility)
      Age <- 10:(n - 1)
      
    }
    
    # 5-year groups
    if(output == "5-year") { 
      
      Age <- seq(from = 10, length.out = n, by = 5)
      
    }
    
    rownames(Asfrmat) <- Age
    
    return(Asfrmat)
    
  }
  
  # -------------------------------
  # CASE 2: download from WPP
  # -------------------------------
  if(is.null(location)) {
    
    stop("You need to provide a location to download the data for Asfrmat")
    
  }
  
  # ------------------------------------------------------ #
  # installed wpp versions
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  # stop if none
  if(length(installed_wpp) == 0) { 
    
    stop("No wpp package installed.")
    
  }
  # find the lates one
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  
  if(latest_wpp < 2022) {
    
    stop("No single ages are availabe in wpp versions earlier that wpp2022.
           Please update the wpp package.")
    
  }
  
  # download mx1dt data from the latest package version
  data("percentASFR1dt", package = latest_wpp)
  data("tfr1",           package = latest_wpp)
  
  # if any date chosen is less then 1950 or more than wpp version + 1
  if(any(AsfrDatesIn < 1950, AsfrDatesIn > (parse_number(latest_wpp) + 1))) {
    cat(paste0(
      "Careful, choosing beyond range 1950-",
      parse_number(latest_wpp)
    ))
    
  }
  # ------------------------------------------------------ #
  # find location code from the provided location
  # if location is misspelled return NA
  
  if(is.numeric(location)) {
    
    location_code <- location
    
  } else {
    
    location_code <- percentASFR1dt %>%
      as_tibble() %>% 
      filter(.data$name %in% location) %>% 
      select("country_code") %>% 
      unique() %>%
      as.numeric()
    
  }
  
  if(verbose) {
    cat(paste0(
      "Downloading ASFR data for ",
      location,
      ", years ",
      paste(AsfrDatesIn, collapse = ", ")
    ),
    sep = "\n")
  }
  
  age <- unique(percentASFR1dt$age)
  
  # -------------------------------
  # Prepare TFR series filtered by country and year
  # -------------------------------
  tfr <- tfr1 %>% 
    as_tibble() %>% 
    pivot_longer(-c("country_code", "name"),
                 names_to  = "year",
                 values_to = "tfr") %>% 
    filter(.data$country_code %in% location_code,
           .data$year < parse_number(latest_wpp) + 1) %>% 
    mutate(year = as.integer(.data$year))
  
  # -------------------------------
  # OUTPUT OPTION: single-year ASFR
  # -------------------------------
  # For single-year ASFR we use percentASFR1dt (pasfr is percent by single-year age)
  # and convert percentages to rates using TFR. After constructing wide matrices
  # we interpolate the columns (years) to requested AsfrDatesIn.
  if(output == "single") { 
    
    out <- percentASFR1dt %>% 
      as_tibble() %>% 
      filter(.data$country_code %in% location_code,
             .data$year <= max(tfr$year)) %>%
      # join TFR so we can convert percentage-of-TFR to ASFR
      left_join(tfr, by = c("country_code", "name", "year")) %>%
      group_by(.data$country_code, .data$name,.data$year) %>% 
      # create asfr
      # Convert pasfr (percent of TFR) to absolute ASFR for each single age.
      # Important: earlier code divided pasfr by sum(pasfr) — that ensures normalization.
      mutate(asfr = (.data$pasfr / sum(.data$pasfr)) * .data$tfr) %>%
      ungroup() %>% 
      select(-c("pasfr", "tfr")) %>% 
      # wide format
      pivot_wider(names_from  = year, 
                  values_from = asfr) %>%
      select(-"age") %>%
      # group rows by country/name in case there are multiple territories (kept for consistency)
      group_nest(.data$country_code, .data$name) %>%
      # interpolate the columns to the requested AsfrDatesIn
      mutate(data = map(
        data,
        ~ interp(
          .x,
          as.numeric(names(.x)),
          as.numeric(AsfrDatesIn),
          extrap = TRUE,
          method = method,
          ...
        ) %>%
          as_tibble()
      )) %>%
      unnest(.data$data) %>%
      mutate(age = age) %>% 
      select(-c("country_code", "name")) %>% 
      column_to_rownames("age") %>% 
      as.matrix()
    
  }
  
  # -------------------------------
  # OUTPUT OPTION: 5-year aggregated ASFR
  # -------------------------------
  # For 5-year groups we aggregate percentASFR into 5-year bins first (by summing
  # normalized percentages), then multiply aggregated percentages by TFR to obtain
  # ASFR for each 5-year group, and then interpolate across years.
  if(output == "5-year") {
    
    # age5: the set of 5-year group labels derived from single-year ages
    age5 <- unique((age %/% 5) * 5)
    
    out <- percentASFR1dt %>% 
      as_tibble() %>% 
      filter(.data$country_code %in% location_code,
             .data$year <= max(tfr$year)) %>%
      group_by(.data$country_code, .data$name, .data$year) %>% 
      # create asfr
      # Convert individual single-year percentages to normalized shares
      mutate(pasfr = (.data$pasfr / sum(.data$pasfr)),
             age   = (.data$age %/% 5) * 5) %>% 
      # sum normalized shares across the 5-year bins
      group_by(.data$country_code, .data$name, .data$year, .data$age) %>%
      summarise(pasfr = sum(pasfr)) %>%
      # attach TFR and convert to ASFR for 5-year groups
      left_join(tfr, by = c("country_code", "name", "year")) %>% 
      mutate(asfr = .data$pasfr * .data$tfr) %>%
      ungroup() %>% 
      select(-c("pasfr", "tfr")) %>% 
      # wide format
      pivot_wider(names_from  = year, 
                  values_from = asfr) %>%
      select(-"age") %>%
      group_nest(.data$country_code, .data$name) %>%
      # interpolate
      # interpolate to requested AsfrDatesIn
      mutate(data = map(
        data,
        ~ interp(
          .x,
          as.numeric(names(.x)),
          as.numeric(AsfrDatesIn),
          extrap = TRUE,
          method = method,
          ...
        ) %>%
          as_tibble()
      )) %>%
      unnest(.data$data) %>%
      mutate(age = age5) %>% 
      select(-c("country_code", "name")) %>% 
      column_to_rownames("age") %>% 
      as.matrix()
    
  } 
  
  return(out)
  
}

# ASFR tests
# When extrapolation range is 5+ years beyond what we have in data
# results in many zeroes (expected?)
# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = NULL,
#   location    = "Argentina",
#   AsfrDatesIn = 1950:2030,
#   method      = "linear",
#   output      = "5-year" 
# )

# power
# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = NULL,
#   location    = "Brazil",
#   AsfrDatesIn = 1950:1960,
#   method      = "power",
#   output      = "single" 
# )

# exponential option is not working here
# returns an error
# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = NULL,
#   location    = "Spain",
#   AsfrDatesIn = 2020,
#   method      = "exponential",
#   output      = "single" 
# )

# with defined ASFR matrix
# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = ASFR_Arg1[, "2020"],
#   location    = "Spain",
#   AsfrDatesIn = 2020,
#   method      = "linear",
#   output      = "single" 
# )

# or in five years
# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = NULL,
#   location    = "Sweden",
#   AsfrDatesIn = 2020,
#   method      = "linear",
#   output      = "5-year" 
# )

# ASFR_Arg1 <- download_Asfr(
#   Asfrmat     = ASFR_Arg1,
#   location    = "Sweden",
#   AsfrDatesIn = 2020,
#   method      = "linear",
#   output      = "5-year" 
# )


