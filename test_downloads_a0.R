purrr::accumulate(0:10, `+`, .dir = "forward")
rev(rev(cumsum(x = 0:10)))
rev(lt_id_L_T(rev(0:10)))

library(DemoTools)
library(tidyverse)

# nLx        = NULL
# location   = "Argentina"
# gender     = "both"
# nLxDatesIn = 1950:2030
# method     = "linear"
# output     = "5-year"
# 
# refDate  <- 1986
# location <- "Brazil"



z <- download_nLx(nLx        = NULL,
                  location   = "Brazil",
                  gender     = "female",
                  nLxDatesIn = nLxDatesIn,
                  method     = "linear",
                  output     = "5-year",
                  radix      = 1)


nLx        = nLxFemale,
location   = location,
# gender   = "female",
nLxDatesIn = floor(nLxDatesIn),


nLx        = NULL
location   = "Brazil"
gender     = "female"
nLxDatesIn = c(1985, 1978)
method     = "linear"
output     = "5-year"
radix      = 1



#' Download or Process nLx Life Table Values
#'
#' This function either returns user-supplied \code{nLx} values in the desired age structure (single, abridged, or 5-year), or downloads and constructs
#' \code{nLx} data from the most recent installed \code{wpp} package for a given location and gender. It supports interpolation and summarization of life table values across different age formats.
#'
#' @param nLx Optional matrix or data frame of user-supplied \code{nLx} values. If provided, the function will return the data formatted according to \code{output}.
#' @param location Character or numeric code specifying the country or region for which to download data. Required if \code{nLx} is not supplied.
#' @param gender Character, one of \code{"male"}, \code{"female"}, or \code{"both"}.
#' @param nLxDatesIn Numeric vector of reference years to download or align life table data.
#' @param method Character specifying interpolation method to use when \code{nLx} data are interpolated across years. Defaults to \code{"linear"}.
#' @param output Character specifying output age structure: \code{"single"}, \code{"abridged"}, or \code{"5-year"} (default).
#' @param radix Numeric. Starting radix for life table construction.
#' @param ... Additional arguments passed to internal interpolation or life table functions.
#'
#' @details
#' The function will look for the latest available \code{wpp} package (e.g. \code{wpp2022}) on the system, extract mortality data (\code{mx1dt}), interpolate missing values for specified years, and calculate \code{nLx} using single-age life tables.
#'
#' When \code{nLx} is provided manually, no download occurs and only the output structure is adjusted according to \code{output}.
#'
#' @return
#' A numeric matrix of \code{nLx} values by age (rows) and year (columns).
#'
#' @importFrom dplyr filter select mutate group_by summarise across rename
#' @importFrom dplyr group_nest ungroup
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map map2
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom stringr str_remove
#' @importFrom readr parse_number
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example using downloaded data
#' nLx_data <- download_nLx(location = "Germany",
#'                          gender = "female",
#'                          nLxDatesIn = c(2000, 2005, 2010),
#'                          output = "abridged")
#'
#' # Example using user-supplied matrix
#' nLx_data <- matrix(runif(20), nrow = 4)
#' download_nLx(nLx = nLx_data,
#'              nLxDatesIn = c(2000, 2005, 2010, 2015),
#'              output = "5-year")
#' }
#'
#' @export

# updated to handle 5-year output ot single year output (age)
download_nLx <- function(nLx      = NULL,
                         location = NULL,
                         gender   = NULL,
                         nLxDatesIn = NULL,
                         method = "linear",
                         output = "5-year",
                         radix = 1,
                         ...) {
  
  verbose <- getOption("basepop_verbose", TRUE)
  
  # if nLx is provided by user just return back the nLx data
  if(!is.null(nLx)) {
    
    nLx           <- as.matrix(nLx)
    colnames(nLx) <- nLxDatesIn
    n             <- nrow(nLx)
    
    if(output == "single") { 
      
      Age <- 0:(n - 1)
      
    }
    
    if(output == "abridged") { 
        
      Age <- c(0, 1, seq(5, (n - 2) * 5, by = 5))
      
    } 
    
    if(output == "5-year") { 
      
        Age <- seq(from = 0, to = n - 1, by = 5)
      
      }
    
    rownames(nLx) <- Age
    
    return(nLx)
    
  }
  
  # until the end of function
  if(is.null(nLx)) {
    
    # if no location return error
    if(is.null(location)) {
      
      stop("You need to provide a location to download the data for nLx")
      
    }
    
    # since we will need to download data anyway, I have changed the
    # function structure a little bit.
    # We first download data, then we ding if ID is provided,
    # then we do filtering and all other checks.
    # I do not think this will change the performance of the function
    # ------------------------------------------------------ #
    # NEW
    # I check which wpp versions are available on the machine
    # if none we stop
    # then I use the latest available for the data grab
    
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
    data("mx1dt", package = latest_wpp)
    
    # if any date chosen is less then 1950 or more than wpp version + 1
    if(any(nLxDatesIn < 1950, nLxDatesIn > (parse_number(latest_wpp) + 1))) {
      
      cat(paste0("Careful, choosing beyond range 1950-", parse_number(latest_wpp), ". "))
      
    }
    
    if(is.numeric(location)) {
      
      location_code <- location
      
    } else {
      
      # find location code from the provided location
      # if location is misspelled return NA
      location_code <- mx1dt %>%
        filter(.data$name %in% location) %>% 
        select("country_code") %>%
        unique() %>%
        as.numeric()
      
    }
    # ------------------------------------------------------ #
    # if we need message print it
    if(verbose) {
      cat(
        paste0(
          "Downloading nLx data for ",
          location,
          ", years ",
          paste(nLxDatesIn, collapse = ", "),
          ", gender ",
          gender
        ),
        sep = "\n"
      )
    }
    
    # standardize input strings to match what we expect
    sex_code <- ifelse(tolower(gender) == "both",
                       "b",
                       ifelse(
                         tolower(gender) == "female",
                         "f",
                         ifelse(tolower(gender) == "male", "m", NA)
                       ))
    
    # sex standardization
    Sex_mortlaws <- ifelse(sex_code == "b", "total", tolower(gender))
    
    stopifnot(`Invalid sex name, please set it to 'both', 'male' or 'female'` = !is.na(sex_code))
    
    # here some data wrangling going on, including
    # calulation of interp()
    # then lt_sinle Lx calculation
    
    out <- mx1dt %>%
      as_tibble() %>%
      filter(.data$country_code %in% location_code,
             .data$year < parse_number(latest_wpp) + 1) %>%
      pivot_longer(-c("country_code", "name", "year", "age"),
                   names_to  = "sex",
                   values_to = "mx") %>%
      mutate(sex = str_remove(.data$sex, "mx"), 
             sex = tolower(.data$sex)) %>%
      filter(.data$sex %in% sex_code) %>%
      pivot_wider(names_from  = year, 
                  values_from = mx) %>%
      select(-"age") %>%
      group_nest(.data$country_code, .data$name, .data$sex) %>%
      # interpolate
      mutate(data = map(
        data,
        ~ interp(
          .x,
          as.numeric(names(.x)),
          as.numeric(nLxDatesIn),
          extrap = TRUE,
          method = method,
          ...
        ) %>%
          as_tibble()
      )) %>% 
      unnest(.data$data) %>% 
      pivot_longer(-c("country_code", "name", "sex"),
                   names_to  = "year",
                   values_to = "mx") %>% 
      group_nest(.data$country_code, .data$name, .data$sex, .data$year) %>%
      # calculate lifetable
      # SRB???
      mutate(data = map2(.x = data, 
                         .y = sex, ~ lt_single_mx(nMx    = .x$mx,
                                                  radix  = radix,
                                                  Sex    = gender,
                                                  sex    = .y,
                                                  a0rule = "ak"))) %>%
      unnest(.data$data) %>%
      select(age = .data$Age, .data$year, .data$nLx) %>% 
      # wide format
      pivot_wider(names_from  = year, 
                  values_from = nLx)
      
    
    
    # what type of age output is desired
    if(output == "abridged") {      
      
      out <- out %>%
        mutate(age = case_when(
          .data$age == 0        ~ 0,
          .data$age %in% c(1:4) ~ 1,
          TRUE                  ~ (.data$age %/% 5) * 5)) %>%                  
        group_by(.data$age) %>%
        summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop") %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    } else if (output == "5-year") { 
      
      out <- out %>%
        mutate(age = (.data$age %/% 5) * 5) %>%                  
        group_by(.data$age) %>%
        summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop") %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
      } else {
      
      out <- out %>%
        column_to_rownames("age") %>% 
        as.matrix()
      
      }
    
    return(out)
    
  }
}

# Asfrmat     = NULL
# location    = "Argentina"
# AsfrDatesIn = 1950:2030
# 
# 
# Asfrmat = NULL
# location = country_code
# AsfrDatesIn = refDate_start:refDate
# method      = "linear"
# output      = "single"



#' Download or Process Age-Specific Fertility Rates (ASFR)
#'
#' This function returns user-supplied ASFR data formatted to a specific age structure or downloads and constructs ASFR values from the most recent available \code{wpp} package for a given location. It supports interpolation across years and aggregation to 5-year or abridged age groups.
#'
#' @param Asfrmat Optional matrix or data frame of user-supplied ASFR values.If provided, the function returns the formatted ASFR data according to \code{output}.
#' @param location Character or numeric code specifying the country or region for which to download data. Required if \code{Asfrmat} is not supplied.
#' @param AsfrDatesIn Numeric vector of reference years for which ASFR data are to be downloaded or aligned.
#' @param method Character specifying interpolation method to use when ASFR data are interpolated across years. Defaults to \code{"linear"}.
#' @param output Character specifying desired output age structure: \code{"single"}, \code{"abridged"}, or \code{"5-year"} (default).
#' @param ... Additional arguments passed to internal interpolation routines.
#'
#' @details
#' The function locates the latest installed \code{wpp} package (e.g. \code{wpp2022}) and extracts percent age-specific fertility rates (\code{percentASFR1dt}) and total fertility rate (\code{tfr1}) data. It then interpolates the data to match requested years (\code{AsfrDatesIn}) and computes ASFR values. If requested, it aggregates them into 5-year age groups using life table exposure (\code{nLx}) values from \code{\link{download_nLx}}.
#'
#' If \code{Asfrmat} is supplied manually, no download occurs and only age-formatting is applied.
#'
#' @return
#' A numeric matrix of ASFR values by age (rows) and year (columns).
#'
#' @importFrom dplyr filter select mutate group_by summarise across rename left_join
#' @importFrom dplyr group_nest ungroup
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column
#' @importFrom stringr str_remove
#' @importFrom readr parse_number
#' @importFrom rlang .data
#'
#' @seealso [download_nLx]
#'
#' @examples
#' \dontrun{
#' # Example using downloaded data
#' asfr_data <- download_Asfr(location = "France",
#'                            AsfrDatesIn = c(2000, 2005, 2010),
#'                            output = "5-year")
#'
#' # Example using user-supplied matrix
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
  
  # if Asfrmat is provided by user just return back the Asfrmat data
  if(!is.null(Asfrmat)) {
    
    Asfrmat           <- as.matrix(Asfrmat)
    colnames(Asfrmat) <- AsfrDatesIn
    n                 <- nrow(Asfrmat)
    
    if(output == "single") { 
      
      Age <- 0:(n-1)
      
    }
    
    if(output == "abridged") { 
      
      Age <- c(0, 1, seq(5, (n - 2) * 5, by = 5))
      
    } 
    
    if(output == "5-year") { 
      
      Age <- seq(from = 0, to = n - 1, by = 5)
      
    }
    
    rownames(Asfrmat) <- Age
    
    return(Asfrmat)
    
  }
  
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
  
  tfr <- tfr1 %>% 
    as_tibble() %>% 
    pivot_longer(-c("country_code", "name"),
                 names_to  = "year",
                 values_to = "tfr") %>% 
    filter(.data$country_code %in% location_code,
           .data$year < parse_number(latest_wpp) + 1) %>% 
    mutate(year = as.integer(.data$year))
  
  out <- percentASFR1dt %>% 
    as_tibble() %>% 
    filter(.data$country_code %in% location_code,
           .data$year <= max(tfr$year)) %>%
    left_join(tfr) %>% 
    # create asfr
    mutate(asfr = (.data$pasfr / 100) * .data$tfr) %>%
    select(-c("pasfr", "tfr")) %>% 
    # wide format
    pivot_wider(names_from  = year, 
                values_from = asfr) %>%
    select(-"age") %>%
    group_nest(.data$country_code, .data$name) %>%
    # interpolate
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
    mutate(age = age)
  
  if(output == "5-year") {
    
    nLx <- download_nLx(nLx = NULL,
                        location   = location,
                        gender     = "female",
                        nLxDatesIn = AsfrDatesIn,
                        method     = "linear",
                        output     = "single",
                        radix      = 100000) %>%
      as.data.frame() %>% 
      rownames_to_column("age") %>%
      pivot_longer(-c("age"),
                   names_to  = "year",
                   values_to = "nLx") %>% 
      mutate(age = as.numeric(.data$age),
             nLx = .data$nLx)
    
    out <- out %>% 
      select(-c("country_code", "name")) %>%
      pivot_longer(-c(age),
                   names_to  = "year",
                   values_to = "asfr") %>%
      left_join(nLx, by = c("age", "year")) %>% 
      filter(!is.na(.data$asfr)) %>%
      mutate(age = (.data$age %/% 5) * 5) %>%
      group_by(.data$year, .data$age) %>%
      summarise(
        asfr = sum(.data$asfr * .data$nLx) / sum(.data$nLx),
        .groups = "drop") %>% 
      pivot_wider(names_from  = year,
                  values_from = asfr) %>% 
      column_to_rownames("age") %>% 
      as.matrix()    
    
  } else { 
    
    out <- out %>% 
      select(-c("country_code", "name")) %>% 
      column_to_rownames("age") %>% 
      as.matrix()
    
  }
  
  return(out)
  
}


SRB        = NULL
location   = "Argentina",
DatesOut = 1950:1952


#' Download or Process Sex Ratio at Birth (SRB)
#'
#' This function either returns user-supplied SRB values for specified years or downloads country-specific SRB data from the most recent installed \code{wpp} package. Missing years are filled with the mean SRB if necessary.
#'
#' @param SRB Optional numeric value of SRB. If provided, this value will be repeated for all \code{DatesOut}.
#' @param location Character or numeric code specifying the country for which SRB data should be downloaded.
#' @param DatesOut Numeric vector of years for which SRB values are required.
#' @param verbose Logical; if \code{TRUE}, prints progress messages. Default is \code{TRUE}.
#'
#' @details
#' The function checks for the latest installed \code{wpp} package (e.g., \code{wpp2022}) and extracts single-year SRB data (\code{sexRatio1}) for the requested country and years. If the country is unavailable in WPP data, a default SRB value of 1.0486 is used. Missing years are interpolated using the mean of available SRB values.
#'
#' @return
#' A named numeric vector of SRB values for each year in \code{DatesOut}.
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example using downloaded SRB for France
#' srb_data <- download_SRB(location = "France",
#'                          DatesOut = c(2000, 2005, 2010))
#'
#' # Example using a user-supplied SRB
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
  if (!is.null(SRB)) {
    
    SRB <- setNames(rep(SRB, length.out = length(DatesOut)), DatesOut)
    
    return(SRB)
    
  }
  
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  if(length(installed_wpp) == 0) stop("No WPP package installed.")
  
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  
  if(latest_wpp < 2022) {
    
    stop("No single years SRB is availabe in wpp versions earlier that wpp2022.
           Please update the wpp package.")
    
  }
  
  data("sexRatio1", package = latest_wpp)
  
  if(is.numeric(location)) {
    
    location_code <- location
    
  } else {
    
    location_code <- sexRatio1 %>%
      as_tibble() %>%
      filter(.data$name %in% location) %>%
      select(.data$country_code) %>%
      unique() %>%
      as.numeric()
  }
  
  if(is.na(location_code)) {
    
    if(verbose) {  
      
      cat(location, "not available in wpp. Using default SRB:", SRB_default, "\n")
      
    }
    
    return(setNames(rep(SRB_default, length(DatesOut)), DatesOut))
  }
  
  if(verbose) {
    cat(paste0("\nDownloading SRB for ", location, 
               " for years ", round(min(DatesOut),1), " to ", round(max(DatesOut),1), "\n"))
  }
  
  dt <- sexRatio1 %>%
    as_tibble() %>%
    filter(.data$country_code == location_code) %>%
    pivot_longer(-c("country_code", "name"), 
                 names_to  = "year", 
                 values_to = "SRB") %>%
    mutate(year = as.numeric(.data$year)) %>%
    filter(.data$year %in% floor(DatesOut))
  
  years_srb <- dt$year
  SRB       <- setNames(dt$SRB, years_srb)
  
  # Fill missing years with mean SRB (if needed) 
  DatesOut_floor <- floor(DatesOut)
  yrs_present    <- DatesOut_floor %in% years_srb
  
  if(any(!yrs_present)) {
    mean_srb <- mean(SRB[as.character(DatesOut_floor[yrs_present])])
    SRB      <- c(SRB, setNames(rep(mean_srb, sum(!yrs_present)), DatesOut_floor[!yrs_present]))
  }
  
  SRB <- SRB[order(as.numeric(names(SRB)))]
  SRB <- setNames(SRB, DatesOut)
  
  return(SRB)
}


# Andreev–Kingkade a0
a0_AK <- function(m0, sex = "f") {
  
  lt_rule_ak_m0_a0(M0 = m0, Sex = sex)
  
}


# Solve for a0 given L0 
solve_a0_from_L0 <- function(L0_target, l0 = 1, sex = "f") {
  
  # Basic sanity check
  if (L0_target >= l0) {
    
    stop("L0_target must be less than l0 (since deaths must occur in the interval).")
    
  }
  
  # Define function whose root gives the correct l1
  f_root <- function(l1) {
    m0 <- (l0 - l1) / L0_target
    a0 <- a0_AK(m0, sex)
    L0_pred <- l1 + a0 * (l0 - l1)
    L0_pred - L0_target
  }
  
  # Numerically solve for l1
  sol <- uniroot(f_root, lower = 0, upper = l0)
  l1  <- sol$root
  
  # Compute implied m0 and a0
  m0 <- (l0 - l1) / L0_target
  a0 <- a0_AK(m0, sex)
  
  # Return everything useful
  return(list(
    a0        = a0,
    m0        = m0,
    l1        = l1,
    L0_target = L0_target,
    sex       = sex
  ))
}


#' Construct Single-Year Base Population Using Reverse Survival
#'
#' This function constructs a single-year age-specific population for both sexes using the reverse survival method. It combines life table information, fertility data, and sex ratio at birth (SRB) to estimate population counts by age, while optionally allowing user-supplied population or fertility data.
#'
#' @param location Optional character name of the country (used for messages only).
#' @param refDate Reference date (numeric) for the population (e.g., census year).
#' @param Age Optional numeric vector of ages. Must match population length if provided.
#' @param country_code Numeric code of the country (required).
#' @param Females_single Optional numeric vector of female population counts by age.
#' @param Males_single Optional numeric vector of male population counts by age.
#' @param nLxFemale Optional female life table \code{nLx} data.
#' @param nLxMale Optional male life table \code{nLx} data (currently unused).
#' @param nLxDatesIn Optional vector of years corresponding to life table data.
#' @param AsfrMat Optional matrix of age-specific fertility rates (ASFR).
#' @param AsfrDatesIn Optional vector of years corresponding to ASFR data.
#' @param SRB Optional sex ratio at birth (males per female) vector or numeric.
#' @param SRBDatesIn Optional vector of years corresponding to SRB data.
#' @param radix Numeric starting radix for life table calculations. Default is 1.
#' @param verbose Logical; if \code{TRUE}, prints progress messages. Default is \code{TRUE}.
#' @param ... Additional arguments passed to internal functions.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Retrieves population, mortality, fertility, and SRB data from the
#'         latest installed \code{wpp} package if not supplied by the user.
#'   \item Performs reverse survival using life table \code{nLx} values to
#'         estimate female exposures.
#'   \item Interpolates or processes fertility rates and calculates births.
#'   \item Applies SRB to distribute births between males and females.
#'   \item Returns the final single-year population by age along with plots
#'         of female and male population distributions.
#' }
#'
#' @return A named list with three elements:
#' \item{result}{A tibble containing single-year population by age, including
#'   male and female counts.}
#' \item{figureF}{A \code{ggplot2} object showing female population by age.}
#' \item{figureM}{A \code{ggplot2} object showing male population by age.}
#'
#' @importFrom dplyr filter select mutate arrange group_by ungroup right_join left_join summarize
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_line
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' basepop <- basepop_single(
#'   country_code = 276,  # Germany
#'   refDate = 2000
#' )
#' basepop$result      # Population tibble
#' basepop$figureF      # Female population plot
#' basepop$figureM      # Male population plot
#' }
#'
#' @export
#' 


basepop_single <- function(location       = NULL,
                           refDate        = NULL,
                           Age            = NULL,
                           country_code   = NULL,
                           Females_single = NULL,
                           Males_single   = NULL,
                           nLxFemale      = NULL,
                           nLxMale        = NULL,
                           nLxDatesIn     = NULL,
                           AsfrMat        = NULL,
                           AsfrDatesIn    = NULL,
                           SRB            = NULL,
                           SRBDatesIn     = NULL,
                           radix          = 1,
                           verbose        = TRUE,
                           ...) {
  
  # options from DemoTools version
  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))
  
  # Jan 1 2000 female pop;
  # note:  1999 means dec 31, 1999, so we treat as jan 1, 2000.
  # we need these dates to filter period
  refDate       <- dec.date(refDate) - 1
  refDate_start <- refDate - 9
  
  # we need minimum of date and country to run a function
  if(is.null(refDate) | is.null(country_code)) {
    
    stop("At least refDate and country_code should be provided for function to work")
    
  }
  
  
  # here we attach the latest wpp package available
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  if(length(installed_wpp) == 0) {
    
    stop("No wpp package installed.")
    
  }
  
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  suppressPackageStartupMessages(library(latest_wpp, character.only = TRUE))
  
  
  if(parse_number(latest_wpp) < 2022) { 
    
    warning("No single ages are availabe in wpp versions earlier that wpp2022. Consider updating the wpp package or change to five year solution.")
    
  }
  
  # if user did not provide population we can calculate it from wpp
  if(is.null(Females_single) | is.null(Males_single)) {
    
    data("popAge1dt", package = latest_wpp)
    
  }
  
  # this data is for future calculations
  data("mx1dt", package = latest_wpp)
  
  # if user provided Population vector we turn it into tibble
  if(!is.null(Females_single)) {
    
    Females_single <- tibble(pop    = Females_single,
                             year   = refDate + 1,
                             age    = row_number() - 1,
                             cohort = year - age - 1)
    
    
  }
  
  # same for Males
  if(!is.null(Males_single)) {
    
    Males_single <- tibble(pop    = Males_single,
                           year   = refDate + 1,
                           age    = row_number() - 1,
                           cohort = year - age - 1)
    
  }
  
  # if user did not provide Population
  if(is.null(Females_single)) {
    
    Females_single <- popAge1dt |>
      filter(country_code == !!country_code,
             year == refDate) |>
      select(year, age, pop = popF) |>
      mutate(year   = year + 1,
             cohort = year - age - 1)
    
  }
  
  # same operation for males
  if(is.null(Males_single)) {
    
    Males_single <- popAge1dt |>
      filter(country_code == !!country_code,
             year == refDate) |>
      select(year, age, pop = popM) |>
      mutate(year   = year + 1,
             cohort = year - age - 1)
    
  }
  
  # Age setup
  if (!is.null(Age)) {
    stopifnot(is_single(Age))
    stopifnot(length(Age) == nrow(Females_single))
    
  } else {
    
    Age <- as.integer(sort(unique(Females_single$age)))
    
  }
  
  
  if(is.null(radix) & !is.null(nLxFemale)) {
    
    radix <- lt_infer_radix_from_1L0(nLxFemale[1, 1])
    
    if(verbose) {
      
      cat(paste0("Setting radix to value of lx: ", radix, 
                 ". Can be overwritten with the `radix` argument"),
          sep = "\n")
      
    }
    
  }
  
  
  if(is.null(radix) & is.null(nLxFemale)) { 
    
    radix <- 1
    cat(paste0("Setting radix to value of lx: ", radix, 
               ". Can be overwritten with the `radix` argument"),
        sep = "\n")
    
    }
  
  
  # nLxFemale = NULL
  # nLxDatesIn = NULL
  # if user did not
  
  # this should be changed!!!!
  if(!is.null(nLxFemale)) {
    
    mxF <- nLxFemale %>% 
      group_by(year) %>% 
      mutate(n  = age2int(.data$age, OAvalue = Inf),
             ax = ifelse(age == 0, 
                         solve_a0_from_L0(L0  = .data$Lxp[.data$age == 0],
                                          l0  = radix,
                                          sex = "f"), 
                         0.5)) %>% 
      unnest(ax) %>% 
      mutate(qx      = (Lxp - lead(Lxp)) / Lxp,
             lx      = lt_id_q_l(nqx = qx, radix = radix),
             dx      = lx * qx,
             Lx      = lx - dx * (1 - ax),
             mx      = lt_id_qa_m(nqx = qx, nax = ax, AgeInt = n),
             cohort  = year - age - 1,
             qx      = ifelse(age == max(age), 1, qx),
             mx      = ifelse(age == max(age), 1 / ax[age == max(age)], mx)) %>% 
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx    = lt_id_q_l(nqx = qx, radix = radix),
        dx    = lx * qx,
        Lx    = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp   = lt_id_q_l(nqx = qx, radix = radix),
        dxp   = lxp * qx,
        Lxp   = lxp - dxp * (1 - ax),
        SxRev = if_else(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
    
  }
  
  if(is.null(nLxFemale)) {
    
    mxF <- mx1dt |>
      filter(country_code == !!country_code,
             between(year, refDate_start, refDate)) |>
      as_tibble() |>
      select(year, age, mx = mxF) |>
      # could try to warp to PC shape here,
      # but uncertain infants. Maybe using
      # an a0 assumption it'd be doable.
      # need cohorts to structure reverse survival
      mutate(
        cohort = year - age - 1,
        age_int = 1,
        ax = if_else(age == 0, lt_rule_1a0_ak(M0 = mx, Sex = "f"), 0.5),
        qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)
      ) |>
      arrange(cohort, age) |>
      group_by(cohort) |>
      mutate(
        lx = lt_id_q_l(nqx = qx, radix = radix),
        dx = lx * qx,
        Lx = lx - dx * (1 - ax),
        SxRev = Lx / lead(Lx, default = 1)
      ) |>
      ungroup() |>
      arrange(year, age) |>
      group_by(year) |>
      mutate(
        lxp = lt_id_q_l(nqx = qx, radix = radix),
        dxp = lxp * qx,
        Lxp = lxp - dxp * (1 - ax),
        SxRev = if_else(year == max(year), Lxp / lead(Lxp), SxRev)
      ) |>
      ungroup() |>
      filter(age < 100) |>
      select(cohort, year, age, SxRev) |>
      arrange(cohort, -age) |>
      group_by(cohort) |>
      # inflation for reverse-surviving
      mutate(inflation_factor = cumprod(SxRev))
  }
  
  expF <- Females_single |>
    select(cohort, pop) |>
    right_join(mxF, by = join_by(cohort)) |>
    mutate(pop_hat = pop * inflation_factor) |>
    select(year, age, pop = pop_hat) |>
    bind_rows(Females_single |> select(year, age, pop)) |>
    filter(between(age, 15, 50)) |>
    arrange(age, year) |>
    mutate(pop_1p1  = lead(pop),
           exposure = (pop + pop_1p1) / 2) |>
    filter(age < 50)
  
  if(!is.null(AsfrMat)) {
    
    AsfrMat <- AsfrMat %>%
      pivot_longer(-c(country_code, name, age),
                   names_to  = "year",
                   values_to = "asfr") %>%
      mutate(year = as.numeric(year)) %>%
      select(-name)
    
    
  }
  
  if(is.null(AsfrMat)) {
    
    AsfrMat <- download_Asfr(Asfrmat     = NULL,
                             location    = country_code,
                             AsfrDatesIn = refDate_start:refDate,
                             method      = "linear",
                             output      = "single") %>%
      as.data.frame() %>% 
      rownames_to_column("age") %>% 
      pivot_longer(-c(age),
                   names_to  = "year",
                   values_to = "asfr") %>%
      mutate(year = as.numeric(year),
             age  = as.numeric(age))
    
  }
  
  
  Bt <- left_join(expF, AsfrMat , by = join_by(year, age)) |>
    filter(year < refDate + 1) |>
    mutate(Bx = asfr * exposure) |>
    group_by(year) |>
    summarize(B = sum(Bx)) |>
    mutate(age = refDate - year) |>
    select(-year)
  
  if(!is.null(SRB)) {
    
    SRB <- tibble(SRB    = SRB,
                  cohort = refDate_start:refDate)
    
  }
  
  if(is.null(SRB)) {
    
    SRB <- download_SRB(SRB      = NULL,
                        location = country_code,
                        DatesOut = refDate_start:refDate,
                        verbose  = TRUE)
    
    SRB <- tibble(SRB    = SRB,
                  cohort = refDate_start:refDate)
    
  }
  
  # final part
  pop_hat <- mx1dt |>
    filter(country_code == !!country_code) |>
    as_tibble() |>
    select(year, age, mx = mxF) |>
    # could try to warp to PC shape here,
    # but uncertain infants. Maybe using
    # an a0 assumption it'd be doable.
    
    # need cohorts to structure reverse survival
    mutate(cohort = year - age - 1,
           age_int = 1,
           ax = if_else(age == 0,
                        lt_rule_1a0_ak(M0 = mx, Sex = "f"),
                        0.5),
           qx = lt_id_ma_q(nMx = mx, nax = ax, AgeInt = age_int)) |>
    filter(between(cohort, refDate_start,refDate),
           between(year, refDate_start, (refDate + 1)),
           age < 10) |>
    arrange(cohort, age) |>
    group_by(cohort) |>
    mutate(lx = lt_id_q_l(nqx = qx, radix = radix),
           dx = lx * qx,
           Lx = lx - dx * (1-ax)) |>
    filter(year == max(year)) |>
    select(age, Lx, cohort) |>
    left_join(Bt, by = join_by(age)) |>
    left_join(SRB, by = join_by(cohort)) |>
    mutate(pop_hat = Lx * B,
           popF = pop_hat * 1 / (1 + SRB),
           popM = pop_hat * SRB / (1 + SRB)) %>%
    ungroup()
  
  # What is Male data is provided???
  # this should be changed!!!!
  # if(!is.null(nLxMale)) {
  #   
  # }
  
  females <- pop_hat |>
    select(age, pop = popF) |>
    ggplot(aes(x = age, y = pop)) +
    geom_line() +
    geom_line(data = Females_single |> filter(age < 10), color = "red")
  
  males <- pop_hat |>
    select(age, pop = popM) |>
    ggplot(aes(x = age, y = pop)) +
    geom_line() +
    geom_line(data = Males_single |> filter(age < 10), color = "red")
  
  # this can be read in from the data
  # we just can take birth from wpp???? we can generate female exposure based on female population
  # reverse survive wpp and take new exposures
  # utils download function 294 interp_co_download_mortailty demotools
  # gets mx then calculates Sx from it.
  # first and last matrix leslie output will be discounted in case it is not cleran dates
  # take the product of the cohort diagonals year age Sx and then make a cohort variable year - age - 1
  # and then arrange cohort age group by cohort summarize Sx = prod(Sx) - survive births forward to the census
  # use Lx downloader and utild for baasepopsingle
  
  return(list(result = pop_hat),
         figureF = females,
         figureM = males
  )
  
}