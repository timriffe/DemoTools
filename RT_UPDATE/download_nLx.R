#' Download or Construct nLx Life Table Values
#'
#' This function either formats user-supplied \code{nLx} data or downloads, interpolates, 
#' and constructs \code{nLx} values from the most recent installed \code{wpp} package. 
#' It supports returning single-year, abridged, or 5-year life table values for a given 
#' location and gender.
#'
#' @param nLx Optional matrix or data frame of user-supplied \code{nLx} values. 
#'   If provided, these will be formatted and returned according to the chosen \code{output} structure.
#' @param location Character or numeric code specifying the country or region for which 
#'   to download data. Required if \code{nLx} is not provided.
#' @param gender Character; one of \code{"male"}, \code{"female"}, or \code{"both"}. 
#'   Determines which mortality rates are used from the \code{wpp} package.
#' @param nLxDatesIn Numeric vector of reference years to download or interpolate life table data for.
#' @param method Character specifying the interpolation method for missing years. 
#'   Defaults to \code{"linear"}. Other methods may be supported by the internal \code{interp()} function.
#' @param output Character specifying the desired output age structure:
#'   \itemize{
#'     \item \code{"single"} – single-year ages (0, 1, 2, …)
#'     \item \code{"abridged"} – standard abridged ages (0, 1, 5, 10, …)
#'     \item \code{"5-year"} – five-year grouped ages (0–4, 5–9, …)
#'   }
#'   Defaults to \code{"5-year"}.
#' @param radix Numeric; radix used for life table construction (usually 1 or 100,000).
#' @param ... Additional arguments passed to internal interpolation or life table functions.
#'
#' @details
#' The function first checks if user-supplied \code{nLx} data is available. 
#' If provided, it reformats the matrix and assigns appropriate age labels 
#' according to \code{output}.  
#'
#' If \code{nLx} is not provided, the function retrieves mortality rate data (\code{mx1dt}) 
#' from the latest installed \code{wpp} package (version ≥ 2022). It then:
#' \enumerate{
#'   \item Identifies the correct country code from the supplied \code{location}.
#'   \item Interpolates single-year mortality rates across specified \code{nLxDatesIn}.
#'   \item Constructs single-age life tables using \code{DemoTools::lt_single_mx()}.
#'   \item Optionally aggregates these results to abridged or 5-year age groups.
#' }
#'
#' This allows consistent generation of life table exposure values (\code{nLx}) 
#' across age structures and time periods, which are essential inputs to 
#' base population estimation and demographic reconstruction.
#'
#' @return
#' A numeric matrix of \code{nLx} values with ages as rows and years as columns.
#'
#' @seealso
#' \code{\link{lt_single_mx}}, \code{\link{lt_single2abridged}}, \code{\link{interp}}.
#'
#' @importFrom dplyr filter select mutate group_by summarise across rename
#' @importFrom dplyr group_nest ungroup
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map map2
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom stringr str_remove
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # Example using downloaded life table data
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
download_nLx <- function(nLx        = NULL,
                         location   = NULL,
                         gender     = NULL,
                         nLxDatesIn = NULL,
                         method     = "linear",
                         output     = "5-year",
                         radix      = 1,
                         ...) {
  
  verbose <- getOption("basepop_verbose", TRUE)
  
  # if nLx is provided by user just return back the nLx data
  # we assume user provides data as vector or as a matrix
  if(!is.null(nLx)) {
    
    nLx           <- as.matrix(nLx)
    colnames(nLx) <- nLxDatesIn
    n             <- nrow(nLx)
    
    
    # this part here deals with the output variation in age
    # user can choose single, abridged or 5 year
    # based on the user choice, the rownames are correspondingly assigned
    # using nrow of nLx
    if(output == "single") { 
      
      Age <- 0:(n - 1)
      
    }
    
    if(output == "abridged") {
      
      Age <- c(0, 1, seq(5, (n - 2) * 5, by = 5))
      
    } 
    
    if(output == "5-year") {
      
      Age <- seq(from = 0, length.out = n, by = 5)
      
    }
    
    rownames(nLx) <- Age
    
    return(nLx)
    
  }
  
  # this part goes until the end of function and deals with missing nLx
  if(is.null(nLx)) {
    
    # At the very least we need a location data. 
    # if no location data is provided we stop with error
    if(is.null(location)) {
      
      stop("You need to provide a location to download the data for nLx")
      
    }
    
    # I check which wpp versions are available on the machine
    # if no wpp is available we stop and offer to use our function to download one
    # otherwise I use the latest available wpp on the machine (in case > 1) for the mortality data 
    
    # installed wpp versions
    installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
    
    # stop if none
    if(length(installed_wpp) == 0) {
      
      stop("No wpp package installed. Please use function install_wpp to download the wpp package.")
      
    }
    
    # find the lates one
    latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
    
    # UPDATE:
    # we basically want wpp to be at leat 2022. We need this for single year data.
    # Since we have a dedicated wpp install function, we agreed that this is ok to do
    if(latest_wpp < 2022) {
      
      stop("No single ages are availabe in wpp versions earlier that wpp2022. Please update the wpp package.")
      
    }
    
    # download mx1dt data from the latest package version
    # this data contain nmx, age, territory, sex, year
    data("mx1dt", package = latest_wpp)
    
    # if any date chosen are less then 1950 or more than wpp version + 1
    if(any(nLxDatesIn < 1950, nLxDatesIn > (parse_number(latest_wpp) + 1))) {
      
      cat(paste0("Careful, choosing beyond range 1950-", parse_number(latest_wpp), ". "))
      
    }
    
    # User can provide location in either number (country_code)
    # OR as a name e.g. Brazil
    # if it is numeric, live as is
    # if it is a name, find the coresponding country number
    if(is.numeric(location)) {
      
      location_code <- location
      
    } else {
      
      # find location code from the provided location if location is misspelled return NA
      location_code <- mx1dt %>%
        filter(.data$name %in% location) %>% 
        select("country_code") %>%
        unique() %>%
        as.numeric()
      
    }
    # ------------------------------------------------------ #
    # if we need message print what we are trying to do exactly
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
    
    # I take the single year mx data, filter country, sex and years up to the latest wpp 
    # then I sue interp for each combination of country and sex for nLxDatesIn (all years used)
    # and thne I calculate the lt_single_mx for same combination
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
      unnest("data") %>% 
      pivot_longer(-c("country_code", "name", "sex"),
                   names_to  = "year",
                   values_to = "mx") %>% 
      group_nest(.data$country_code, .data$name, .data$sex, .data$year) %>%
      # calculate lifetable
      # SRB???
      mutate(data = map2(.x = data, 
                         .y = sex, ~ lt_single_mx(nMx    = .x$mx,
                                                  radix  = radix,
                                                  Sex    = .y,
                                                  a0rule = "ak")))
    
    # now if the output is single ages, we do not need no extra
    if(output == "single") {
      
      out <- out %>% 
        unnest("data") %>%
        select(age = "Age", "year", "nLx") %>% 
        # wide format
        pivot_wider(names_from  = year, 
                    values_from = nLx) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    }
    
    # if output is abridged, then we calculate abridged LT from single with DemoTools
    if(output == "abridged") {
      
      out <- out %>% 
        mutate(data = map2(.x = data,
                           .y = sex, ~ .x %>%
                             reframe(lt_single2abridged(lx  = lx,
                                                        ex  = ex,
                                                        nLx = nLx)) %>% 
                             select("Age", "nLx"))) %>%
        unnest("data") %>% 
        select(age = "Age", "year", "nLx") %>% 
        # wide format
        pivot_wider(names_from  = year, 
                    values_from = nLx) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    } 
    
    # finally if output is 5-years we use abridged LT, but we add first two ages  
    if(output == "5-year") { 
      
      out <- out %>% 
        mutate(data = map2(.x = data,
                           .y = sex, ~ .x %>%
                             reframe(lt_single2abridged(lx  = lx,
                                                        ex  = ex,
                                                        nLx = nLx)) %>% 
                             select("Age", "nLx"))) %>%
        unnest("data") %>% 
        select(age = "Age", "year", "nLx") %>% 
        # wide format
        pivot_wider(names_from  = year, 
                    values_from = nLx) %>%
        mutate(age = ifelse(age %in% c(0, 1), 0, age)) %>% 
        group_by(.data$age) %>%
        summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .groups = "drop") %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    }
    
    # Done
    return(out)
    
  }
}

library(tidyverse)
library(DemoTools)

# Here I test download_nLx function with different arguments
# lets try variate sex, country and nLxDatesIn and radix adn method
# all good, all works
# Lxb <- download_nLx(
#   location   = "Spain",
#   gender     = "both",
#   nLxDatesIn = c(1980, 1985),
#   method     = "linear",
#   output     = "abridged",
#   radix      = 1
# )

# Lxf <- download_nLx(
#   location   = "Brazil",
#   gender     = "female",
#   nLxDatesIn = c(1990:2000),
#   method     = "exponential",
#   output     = "abridged",
#   radix      = 100
# )

# power works strangely
# Lxm <- download_nLx(
#   location   = "Argentina",
#   gender     = "male",
#   nLxDatesIn = c(2000:2022),
#   method     = "power",
#   output     = "abridged",
#   radix      = 100000
# )

# Now lets variate the output. 
# Single
# All works
# LxmSingle <- download_nLx(
#   location   = "Argentina",
#   gender     = "male",
#   nLxDatesIn = c(2000:2022),
#   method     = "linear",
#   output     = "single",
#   radix      = 100000
# )

# Abridged
# LxmAbridged <- download_nLx(
#   location   = "Argentina",
#   gender     = "male",
#   nLxDatesIn = c(2000:2022),
#   method     = "linear",
#   output     = "abridged",
#   radix      = 100000
# )

# 5-year
# LxmFive_year <- download_nLx(
#   nLx        = NULL,
#   location   = "Argentina",
#   gender     = "male",
#   nLxDatesIn = c(2000:2022),
#   method     = "linear",
#   output     = "5-year",
#   radix      = 100000
# )

# looks ok to me
# plot(LxmSingle[,"2021"])
# plot(LxmAbridged[,"2021"])
# plot(LxmFive_year[,"2021"])

# finally. lets try what happens when user provides nLx himself
# IMPORTANT NOTE: if user provides data in single years, but asks for output in 5 or abridged
# it will return the DT but rownames will be strange
# download_nLx(nLx        = LxmSingle[,"2021"],
#              output     = "single",
#              nLxDatesIn = 2021)
# download_nLx(nLx = LxmAbridged[,"2021"],
#              output     = "abridged",
#              nLxDatesIn = 2021)
# download_nLx(nLx = LxmFive_year[,"2021"],
#              output     = "5-year",
#              nLxDatesIn = 2021)

# will return bad rownames (age)
# download_nLx(nLx    = LxmFive_year[,"2021"],
#              output = "single")
