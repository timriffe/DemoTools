
# These utils might be used by basepop, interp_coh, OPAG, mig_resid*,
# and potentially others.

# Authors: Tim, Rustam
#' Download or Construct nLx Life Table Values
#'
#' This function either formats user-supplied \code{nLx} data or downloads, interpolates, and constructs \code{nLx} values from the most recent installed \code{wpp} package. It supports returning single-year, abridged, or 5-year life table values for a given location and gender.
#'
#' @param nLx Optional matrix or data frame of user-supplied \code{nLx} values. If provided, these will be formatted and returned according to the chosen \code{output} structure.
#' @param location Character or numeric code specifying the country or region for which to download data. Required if \code{nLx} is not provided.
#' @param gender Character; one of \code{"m"}, \code{"f"}. Determines which mortality rates are used from the \code{wpp} package.
#' @param nLxDatesIn Numeric vector of reference years to download or interpolate life table data for.
#' @param Age Numeric vector of ages Age that are present in `nLx`
#' @param method Character specifying the interpolation method for missing years. Defaults to \code{"linear"}. Other methods may be supported by the internal \code{interp()} function.
#' @param output Character specifying the desired output age structure:
#'   \itemize{
#'     \item \code{"single"}: single-year ages (0, 1, 2, …)
#'     \item \code{"abridged"}: standard abridged ages (0, 1, 5, 10, …)
#'     \item \code{"5-year"}: five-year grouped ages (0–4, 5–9, …)
#'   }
#'   Defaults to \code{"5-year"}.
#' @param radix Numeric; radix used for life table construction (usually 1 or 100,000), defaults to 1.
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
#' \code{\link{lt_single_mx}}, 
#' \code{\link{lt_single2abridged}}, 
#' \code{\link{interp}}.
#' \code{\link{graduate_pclm}}.
#'
#' @importFrom dplyr select mutate summarise across rename group_nest full_join
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map map2
#' @importFrom tibble column_to_rownames as_tibble
#' @importFrom stringr str_remove
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom utils data
#' @importFrom tidyselect where
#'
#' @examples
#' # Example using downloaded life table data
#' nLx_data <- downloadnLx(location = "Germany",
#'                          gender = "female",
#'                          nLxDatesIn = c(2000, 2005, 2010),
#'                          output = "abridged")
#'
#' # Example using user-supplied matrix
#' nLx_data <- matrix(runif(20), nrow = 4)
#' downloadnLx(nLx = nLx_data,
#'              nLxDatesIn = c(2000, 2005, 2010, 2015, 2020),
#'              output = "5-year")
#'
#' @export

# updated to handle 5-year output ot single year output (age)
downloadnLx <- function(nLx        = NULL,
                        location   = NULL,
                        gender     = NULL,
                        nLxDatesIn = NULL,
                        Age        = NULL,
                        method     = "linear",
                        output     = "5-year",
                        radix      = 1,
                        ...) {
  
  if(is.null(Age)) {
    
  prt_msg <- FALSE
  
  } else { 
    
    prt_msg <- TRUE
    
    }
  output <- match.arg(output, c("single", "abridged", "5-year"))
  
  verbose <- getOption("basepop_verbose", TRUE)
  
  # If nLx is given and Age is also given then nrow() has to equal length(). 
  if(!is.null(nLx) && !is.null(Age) && nrow(nLx) != length(Age)) {
    stop(
      "Inconsistent input: when both `nLx` and `Age` are provided, ",
      "`nrow(nLx)` must equal `length(Age)`. \n",
      "Currently: nrow(nLx) = ", nrow(nLx),
      ", length(Age) = ", length(Age), "."
    )
  }
  
  # check if nLxDatesIn are provided and nLx too. Cols of nLx should match nLxDatesIn
  if(!is.null(nLx) && !is.null(nLxDatesIn) && ncol(nLx) != length(nLxDatesIn)) {
    stop(
      "Inconsistent input: when both `nLx` and `nLxDatesIn` are provided, ",
      "`ncol(nLx)` must equal `length(nLxDatesIn)`. \n",
      "Currently: ncol(nLx) = ", ncol(nLx),
      ", length(nLxDatesIn) = ", length(nLxDatesIn), "."
    )
  }
  
  
  if(!is.null(nLx)) {
    
    nLx <- as.matrix(nLx)
    
    # If nLxDatesIn is provided we reassign column names
    if(!is.null(nLxDatesIn)) {
      
      colnames(nLx) <- nLxDatesIn
      
    }
    
    # Now if Age is not provided we assume age with names2age() function
    if(is.null(Age)) {
      
      Age <- names2age(rownames(nLx))
      
      # if age is provided we reassign rownames names
    } else if(!is.null(Age)) {
      
      rownames(nLx) <- Age
      
    }
    
    return(nLx)
  }
  
  
  # this part goes until the end of function and deals with missing nLx
  if(is.null(nLx)) {
    
    # At the very least we need a location data. 
    # if no location data is provided we stop with error
    if(is.null(location)) {
      
      stop("You need to provide a location to download the data for nLx. Either UN numeric code or official UN location name. You can check the names and codes here https://en.wikipedia.org/wiki/ISO_3166-1_numeric. NOTE: If the code starts with zeroes they should be avoided such as 004, should be read as 4.")
      
    }
    
    # I check which wpp versions are available on the machine
    # if no wpp is available we stop and offer to use our function to download one
    # otherwise I use the latest available wpp on the machine (in case > 1) for the mortality data 
    
    # installed wpp versions
    installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
    
    # stop if none
    if(length(installed_wpp) == 0) {
      
      stop("No wpp package installed. Please use the function chack_and_load_latest_wpp() to install a wpp package version.")
      
    }
    
    # find the lates one
    latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
    
    
    # if any date chosen are less then 1950 or more than wpp version + 1
    if(any(nLxDatesIn < 1950, nLxDatesIn > (parse_number(latest_wpp) + 1))) {
      
      message(paste0("Careful, choosing beyond range 1950-", parse_number(latest_wpp), ". "))
      
    }
    
    # ------------------------------------------------------ #
    # if we need message print what we are trying to do exactly
    if(verbose) {
      message(
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
    sex_code <- ifelse(tolower(gender)== "male",   "m", gender)
    sex_code <- ifelse(tolower(gender)== "men",    "m", sex_code)
    sex_code <- ifelse(tolower(gender)== "female", "f", sex_code)
    sex_code <- ifelse(tolower(gender)== "women",  "f", sex_code)
    sex_code <- ifelse(tolower(gender)== "w",      "f", sex_code)
    sex_code <- ifelse(tolower(sex_code)%in% c("m", "f"),  sex_code, NA_character_)
    
    if(is.na(sex_code)) { 
      
      stop("Invalid sex input. Please use male (men), or female (women), as an input, or use only first letter i.e. m, f (w).")
      
    }    
    
    # UPDATE:
    # we basically want wpp to be at leat 2022. We need this for single year data.
    # Since we have a dedicated wpp install function, we agreed that this is ok to do
    if(parse_number(latest_wpp) < 2022) {
      
      warning("No single ages are available in wpp versions earlier that wpp2022. Please update the wpp package using check_and_load_latest_wpp function.")
      
      env <- new.env()      
      # in this case we need mxM, mxF, popM, and popF
      data("mxM",  package = latest_wpp, envir = env)
      data("mxF",  package = latest_wpp, envir = env)
      # population per 1000 
      data("popF", package = latest_wpp, envir = env)
      data("popM", package = latest_wpp, envir = env)
      
      mxM  <- env$mxM
      mxF  <- env$mxF
      popF <- env$popF
      popM <- env$popM
      
      # define location
      if(is.numeric(location)) {
        
        location_code <- location
        
      } else {
        
        # find location code from the provided location if location is misspelled return NA
        location_code <- mxM %>%
          dplyr::filter(.data$name %in% location) %>% 
          select("country_code") %>%
          unique() %>%
          as.numeric()
        
      }
      
      # years for  interpolation
      anchor_years <- get_bracketing_years(Dates = nLxDatesIn)
      
      # prepare population data
      popM <- popM %>% 
        dplyr::filter(.data$country_code %in% location_code) %>% 
        mutate("sex" = "m") %>% 
        pivot_longer(-c("country_code", "name", "sex", "age"),
                     names_to  = "year",
                     values_to = "pop") %>% 
        mutate("year" = parse_number(.data$year)) %>%
        dplyr::filter(.data$year %in% anchor_years)
      
      popF <- popF %>% 
        dplyr::filter(.data$country_code %in% location_code) %>% 
        mutate("sex" = "f") %>% 
        pivot_longer(-c("country_code", "name", "sex", "age"),
                     names_to  = "year",
                     values_to = "pop") %>% 
        mutate("year" = parse_number(.data$year)) %>%
        dplyr::filter(.data$year %in% anchor_years)
      
      # pop
      population <- popM %>% 
        full_join(popF, by = c("country_code", "name", "age",
                               "sex", "year", "pop")) %>% 
        mutate("age" = parse_number(.data$age),
               # population is per 1000
               "pop" = .data$pop * 1000)
      
      # prepare mx data
      mxM <- mxM %>% 
        dplyr::filter(.data$country_code %in% location_code) %>%
        mutate("sex" = "m") %>% 
        pivot_longer(-c("country_code", "name", "sex", "age"),
                     names_to  = "year",
                     values_to = "mx") %>% 
        mutate("year" = parse_number(.data$year)) %>%
        dplyr::filter(.data$year %in% anchor_years)
      
      mxF <- mxF %>% 
        dplyr::filter(.data$country_code %in% location_code) %>%
        mutate("sex" = "f") %>% 
        pivot_longer(-c("country_code", "name", "sex", "age"),
                     names_to  = "year",
                     values_to = "mx") %>% 
        mutate("year" = parse_number(.data$year)) %>%
        dplyr::filter(.data$year %in% anchor_years)
      
      # A problem the ages in mx are 0, 1-4, ...
      # in Population they are 0-4, ... so we need to combine mortality me thinks
      # the safest way is through the lifetable me thinks
      # first prepare mx data
      mx <- mxM %>% 
        full_join(mxF, by = c("country_code", "name", "age",
                              "sex", "year", 'mx')) %>%
        group_nest(.data$country_code, .data$name, .data$sex, .data$year) %>% 
        mutate(data = map2(.x = .data$data, 
                           .y = .data$sex, ~ lt_abridged(nMx    = .x$mx,
                                                         radix  = radix,
                                                         Sex    = .y,
                                                         Age    = .x$age,
                                                         a0rule = "ak") %>% 
                             select("Age", "ndx", "nLx") %>%
                             mutate(Age = base::replace(.data$Age, .data$Age %in% 0:1, 0)) %>%
                             summarise(
                               ndx = sum(.data$ndx),
                               nLx = sum(.data$nLx),
                               nMx = .data$ndx / .data$nLx,
                               .by = "Age"
                             ) %>% 
                             select(-c("ndx", "nLx")))) %>% 
        unnest("data")
      
      # now that age groups match we can graduate
      mx1 <- population %>% 
        full_join(mx, by = c("country_code", "name", "age" = "Age", "sex", "year")) %>% 
        mutate("deaths" = .data$nMx * .data$pop) %>% 
        group_nest(.data$country_code, .data$name, .data$year, .data$sex) %>% 
        mutate(data1 = map(data, ~ graduate_pclm(Value  = .x$deaths,
                                                 Age    = .x$age,
                                                 OAG    = TRUE,
                                                 offset = .x$pop)%>% 
                             tibble(age = 0:(length(.) - 1), 
                                    mx = .))) %>% 
        select(-"data") %>% 
        unnest("data1")
      
      out <- mx1 %>%
        # Here filter sex
        dplyr::filter(.data$sex %in% sex_code) %>%
        pivot_wider(names_from  = "year", 
                    values_from = "mx") %>%
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
        mutate(data = map2(.x = .data$data, 
                           .y = .data$sex, ~ lt_single_mx(nMx    = .x$mx,
                                                          radix  = radix,
                                                          Sex    = .y,
                                                          a0rule = "ak")))
      
      if(output == "single") {
        
        out <- out %>% 
          unnest("data") %>%
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>%
          column_to_rownames("age") %>% 
          as.matrix()
        
      }
      
      # if output is abridged, then we calculate abridged LT from single with DemoTools
      if(output == "abridged") {
        
        out <- out %>% 
          mutate(data = map2(.x = .data$data,
                             .y = .data$sex, ~ .x %>% 
                               reframe(lt_single2abridged(lx  = lx,
                                                          ex  = ex,
                                                          nLx = nLx)) %>% 
                               select("Age", "nLx"))) %>%
          unnest("data") %>% 
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>%
          column_to_rownames("age") %>% 
          as.matrix()
        
      } 
      
      # finally if output is 5-years we use abridged LT, but we add first two ages  
      if(output == "5-year") { 
        
        out <- out %>% 
          mutate(data = map2(.x = .data$data,
                             .y = .data$sex, ~ .x %>% 
                               reframe(lt_single2abridged(lx  = lx,
                                                          ex  = ex,
                                                          nLx = nLx)) %>% 
                               select("Age", "nLx"))) %>%
          unnest("data") %>% 
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>%
          mutate("age" = replace(.data$age, .data$age %in% 0:1, 0)) %>%
          summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), 
                    .by = "age") %>% 
          column_to_rownames("age") %>% 
          as.matrix()
        
      }
      
      # now if wpp version is 2022+ !!!!
    } else {
      
      env <- new.env()
      # download mx1dt data from the latest package version
      # this data contain nmx, age, territory, sex, year
      data("mx1dt", package = latest_wpp, envir = env)
      mx1dt <- env$mx1dt
      
      # User can provide location in either number (country_code)
      # OR as a name e.g. Brazil
      # if it is numeric, live as is
      # if it is a name, find the coresponding country number
      if(is.numeric(location)) {
        
        location_code <- location
        
      } else {
        
        # find location code from the provided location if location is misspelled return NA
        location_code <- mx1dt %>%
          dplyr::filter(.data$name %in% location) %>% 
          select("country_code") %>%
          unique() %>%
          as.numeric()
        
      }
      
      # I take the single year mx data, filter country, sex and years up to the latest wpp 
      # then I sue interp for each combination of country and sex for nLxDatesIn (all years used)
      # and thne I calculate the lt_single_mx for same combination
      out <- mx1dt %>%
        dplyr::filter(.data$country_code %in% location_code,
                      .data$year < parse_number(latest_wpp) + 1) %>%
        pivot_longer(-c("country_code", "name", "year", "age"),
                     names_to  = "sex",
                     values_to = "mx") %>%
        mutate("sex" = str_remove(.data$sex, "mx"), 
               "sex" = tolower(.data$sex)) %>%
        dplyr::filter(.data$sex %in% sex_code) %>%
        pivot_wider(names_from  = "year", 
                    values_from = "mx") %>%
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
        mutate(data = map2(.x = .data$data, 
                           .y = .data$sex, ~ lt_single_mx(nMx    = .x$mx,
                                                          radix  = radix,
                                                          Sex    = .y,
                                                          a0rule = "ak")))
      
      # now if the output is single ages, we do not need no extra
      if(output == "single") {
        
        out <- out %>% 
          unnest("data") %>%
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>% 
          column_to_rownames("age") %>% 
          as.matrix()
        
      }
      
      # if output is abridged, then we calculate abridged LT from single with DemoTools
      if(output == "abridged") {
        
        out <- out %>% 
          mutate(data = map2(.x = .data$data,
                             .y = .data$sex, ~ .x %>% 
                               reframe(lt_single2abridged(lx  = lx,
                                                          ex  = ex,
                                                          nLx = nLx)) %>% 
                               select("Age", "nLx"))) %>%
          unnest("data") %>% 
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>% 
          column_to_rownames("age") %>% 
          as.matrix()
        
      } 
      
      # finally if output is 5-years we use abridged LT, but we add first two ages  
      if(output == "5-year") { 
        
        out <- out %>% 
          mutate(data = map2(.x = .data$data,
                             .y = .data$sex, ~ .x %>% 
                               reframe(lt_single2abridged(lx  = lx,
                                                          ex  = ex,
                                                          nLx = nLx)) %>% 
                               select("Age", "nLx"))) %>%
          unnest("data") %>% 
          select("age" = "Age", "year", "nLx") %>% 
          # wide format
          pivot_wider(names_from  = "year", 
                      values_from = "nLx") %>%
          mutate("age" = ifelse(.data$age %in% c(0, 1), 0, .data$age)) %>% 
          summarise(across(where(is.numeric), ~ sum(., na.rm = TRUE)), .by = "age") %>% 
          column_to_rownames("age") %>% 
          as.matrix()
        
      }
      
    }
    
    # Check and Assign user provided ages if they comply with requirements
    if(prt_msg & length(Age) != nrow(out)) { 
      
      warning("The ages you provided do not match with the requested data. Assigning ages automatically from the wpp data")
      
    } else { 
      
      message("Using the provided Age vector for rownames of resulting nLx matrix")
      
    }
    
    return(out) 
    
  }
}


# Extract ASFR estimates from WPP2019. Mainly an util function for other ones.
# @description We extract `ASFRx` from `wpp2019`, interpolated to exact dates. Different methods available.
# A vector of countries can handle, but with an unique sex. Row names are not indicative of countries.
# @param Asfrmat numeric.
# @param location vector. UN Pop Div `LocName` or `LocID`
# @param AsfrDatesIn numeric. Vector of decimal dates.
# @param method character. Could be `"linear"`, `"exponential"`, or `"power"`
# 
# @return numeric matrix interpolated asfr
# @export
# @importFrom fertestr get_location_code
# @importFrom fertestr is_LocID
# @importFrom stats setNames
# @examples
# # Total fertility ratio calculated from ASFRx downloaded from WPP19.
# # See `downloadnLx` for analogous examples on multiple countries or using codes instead of names.
# ASFR_Arg <- downloadAsfr(Asfrmat = NULL, location = "Argentina", AsfrDatesIn = 1950:2025)
# \dontrun{
# plot(1950:2025, as.numeric(colSums(ASFR_Arg))*5, xlab = "Year", ylab="TFR", ylim=c(1.5,4), t="l")
# }
# downloadAsfr <- function(Asfrmat, location = NULL, AsfrDatesIn, method="linear") {
# 
#   verbose <- getOption("basepop_verbose", TRUE)
# 
#   if (!is.null(Asfrmat)) {
#     # TR: can we assume colnames are AsfrDatesIn ?
#     return(Asfrmat)
#   }
# 
#   # stop/warnings
#   if (is.null(location)){
#     stop("You need to provide a location to download the data for Asfrmat")
#   }
#   if (!any(fertestr::is_LocID(location))) {
#     location_code <- fertestr::get_location_code(location)
#   }else {
#     location_code <- as.integer(location)
#   }
#   if (verbose) {
#     cat(paste0("Downloading ASFR data for ", location, ", years ", paste(AsfrDatesIn,collapse=", ")), sep = "\n")
#   }
#   if(any(AsfrDatesIn<1950,AsfrDatesIn>2025)){
#     cat("Careful, extrapolating beyond range 1950-2025")
#   }
# 
#   # initial data
#   asfr_wpp19    <-DemoToolsData::WPP2019_asfr
# 
#   # spread format
#   asfr_ctry     <- asfr_wpp19[asfr_wpp19$LocID %in% location_code,] %>%
#                       as.data.frame() %>%
#                       stats::reshape(direction = "wide", idvar = c("LocID","AgeStart"),
#                               timevar = "Year", v.names = "ASFR")
# 
#   # interp/extrap
#   out <- interp(asfr_ctry[,-c(1:3)], seq(1953,2023,5),
#                 as.numeric(AsfrDatesIn),
#                 extrap = TRUE, method = method) %>%
#                 as.data.frame() %>%
#                 stats::setNames(as.character(AsfrDatesIn)) %>%
#                 as.matrix()
# 
#   # combination as rowname
#   rownames(out) <- asfr_ctry$AgeStart
# 
#   return(out)
# }

# Authors: Tim, Rustam
#' Download or Process Age-Specific Fertility Rates (ASFR)
#'
#' This function either returns user-supplied ASFR data formatted to a specific age
#' structure or downloads and constructs ASFR values from the most recent installed
#' \code{wpp} package for a given location. It supports interpolation across years
#' and aggregation to 5-year age groups.
#'
#' When \code{Asfrmat} is supplied, the function:
#' \enumerate{
#'   \item Converts it to a matrix (if needed),
#'   \item Assigns row names (ages) and column names (years),
#'   \item Returns the formatted ASFR matrix.
#' }
#'
#' When \code{Asfrmat} is not supplied, the function:
#' \enumerate{
#'   \item Finds the latest installed \code{wpp} package (>= 2022 for single-year data),
#'   \item Loads \code{percentASFR1dt} and \code{tfr1} datasets,
#'   \item Constructs ASFR for single-year or 5-year age groups,
#'   \item Interpolates ASFR across requested years \code{AsfrDatesIn} using \code{interp()}.
#' }
#'
#' @param Asfrmat Optional numeric matrix or data frame of ASFR values (rows = ages, columns = years).
#' @param location Character or numeric code specifying the country/region to download ASFR for.
#' @param AsfrDatesIn Numeric vector of years to interpolate or extract ASFR for.
#' @param Age Optional numeric vector of ages corresponding to rows of \code{Asfrmat}.
#' @param method Character specifying interpolation method for \code{interp()} (default = "linear").
#' @param output Character specifying age structure of output: \code{"single"} (single-year) or \code{"5-year"} (default).
#' @param ... Additional arguments passed to the interpolation function \code{interp()}.
#'
#' @return A numeric matrix of ASFR values:
#' \itemize{
#'   \item Rows: age groups (single-year or 5-year bins)
#'   \item Columns: years (\code{AsfrDatesIn})
#' }
#'
#' @details
#' For single-year output, ASFR is calculated by multiplying the percentage ASFR
#' from WPP by TFR and interpolating to requested years. For 5-year output, percentages
#' are aggregated into 5-year bins before multiplying by TFR and interpolating.
#' If \code{Asfrmat} is supplied by the user, it is returned after formatting.
#'
#' @seealso 
#' \code{\link{interp}}, 
#' \code{\link{graduate_pclm}}
#' @importFrom dplyr select mutate summarise left_join group_nest rename reframe full_join
#' @importFrom tidyr pivot_longer pivot_wider unnest
#' @importFrom purrr map
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom utils data 
#' 
#' @examples
#' # Using user-supplied ASFR matrix
#' ASFR_user <- matrix(runif(5*3), nrow = 5, ncol = 3)
#' rownames(ASFR_user) <- c(15, 20, 25, 30, 35)
#' colnames(ASFR_user) <- c(2000, 2005, 2010)
#' ASFR_mat <-  downloadAsfr(Asfrmat = ASFR_user)
#'
#' # Download single-year ASFR for Argentina
#' ASFR_Arg <- downloadAsfr(
#'   Asfrmat = NULL,
#'   location = "Argentina",
#'   AsfrDatesIn = 1950:1960,
#'   method = "linear",
#'   output = "single"
#' )
#'
#'ASFR_Arg1 <- downloadAsfr(
#'  Asfrmat     = NULL,
#'  location    = "Argentina",
#'  AsfrDatesIn = 1950:1960,
#'  method      = "linear",
#'  output      = "single"
#')
#'
#'
#'ASFR_Arg2 <- downloadAsfr(
#'  Asfrmat     = NULL,
#'  location    = "Argentina",
#'  AsfrDatesIn = 1950:1963,
#'  method      = "linear",
#'  output      = "5-year"
#')
#' @export

downloadAsfr <- function(Asfrmat     = NULL,
                         location    = NULL,
                         AsfrDatesIn = NULL,
                         Age         = NULL,
                         method      = "linear",
                         output      = "5-year",
                         ...) {
  
  if(is.null(Age)) {
    
    prt_msg <- FALSE
    
  } else { 
    
    prt_msg <- TRUE
    
  }
  
  verbose <- getOption("basepop_verbose", TRUE)
  
  # If Asfrmat is given and Age is also given then nrow() has to equal length(). 
  if(!is.null(Asfrmat) && !is.null(Age) && nrow(Asfrmat) != length(Age)) {
    stop(
      "Inconsistent input: when both `Asfrmat` and `Age` are provided, ",
      "`nrow(Asfrmat)` must equal `length(Age)`. \n",
      "Currently: nrow(Asfrmat) = ", nrow(Asfrmat),
      ", length(Age) = ", length(Age), "."
    )
  }
  
  # check if AsfrDatesIn are provided and Asfrmat too. Cols of Asfrmat should match AsfrDatesIn
  if(!is.null(Asfrmat) && !is.null(AsfrDatesIn) && ncol(Asfrmat) != length(AsfrDatesIn)) {
    stop(
      "Inconsistent input: when both `Asfrmat` and `AsfrDatesIn` are provided, ",
      "`ncol(Asfrmat)` must equal `length(AsfrDatesIn)`. \n",
      "Currently: ncol(Asfrmat) = ", ncol(Asfrmat),
      ", length(AsfrDatesIn) = ", length(AsfrDatesIn), "."
    )
  }
  
  # -------------------------------
  # CASE 1: user-supplied Asfrmat
  # -------------------------------
  # If the user provides Asfrmat we simply coerce to matrix, assign column names
  # (years) and set rownames (ages) according to requested output format.
  
  if(!is.null(Asfrmat)) {
    
    Asfrmat <- as.matrix(Asfrmat)
    
    # If AsfrDatesIn is provided we reassign column names
    if(!is.null(AsfrDatesIn)) {
      
      colnames(Asfrmat) <- AsfrDatesIn
      
    }
    
    # Now if Age is not provided we assume age with names2age() function
    if(is.null(Age)) {
      
      Age <- names2age(rownames(Asfrmat))
      
      # if age is provided we reassign rownames names
    } else if(!is.null(Age)) {
      
      rownames(Asfrmat) <- Age
      
    }
    
    return(Asfrmat)
  }
  
  # -------------------------------
  # CASE 2: download from WPP
  # -------------------------------
  # if location is not provided stop and show where and how to retrive it
  if(is.null(location)) {
    
    stop("You need to provide a location to download the data for Asfrmat. Either UN numeric code or official UN location name. You can check the names and codes here https://en.wikipedia.org/wiki/ISO_3166-1_numeric. NOTE: If the code starts with zeroes they should be avoided such as 004, should be read as 4.")
    
  }
  
  # ------------------------------------------------------ #
  # installed wpp versions
  installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
  
  # stop if none
  if(length(installed_wpp) == 0) { 
    
    stop("No wpp package installed. Please use the function chack_and_load_latest_wpp() to install a wpp package version.")
    
  }
  # find the lates one
  latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
  
  # if any date chosen is less then 1950 or more than wpp version + 1
  if(any(AsfrDatesIn < 1950, AsfrDatesIn > (parse_number(latest_wpp) + 1))) {
    message(paste0(
      "Careful, choosing beyond range 1950-",
      parse_number(latest_wpp)
    ))
    
  }
  
  if(verbose) {
    message(paste0(
      "Downloading ASFR data for ",
      location,
      ", years ",
      paste(AsfrDatesIn, collapse = ", ")
    ),
    sep = "\n")
  }
  
  # in case single user requested single data output from wpp < 2022 we graduate
  if(parse_number(latest_wpp) < 2022) {
    
    warning("No single ages or years are availabe in wpp versions earlier than wpp2022.
            Please update the wpp package. Currently the pclm-graduated data will be provided.")
    
    # create environment for wpp dta
    env <- new.env()    
    
    # in this case we need tfr, percentASFR, and popF
    data("percentASFR", package = latest_wpp, envir = env)
    data("tfr",         package = latest_wpp, envir = env)
    data("popF",        package = latest_wpp, envir = env)
    
    percentASFR <- env$percentASFR
    tfr         <- env$tfr
    popF        <- env$popF
    
    # define location. if it is numeric use as is. 
    # if it is character say "USA" find the corresponding numeric code from data
    if(is.numeric(location)) {
      
      location_code <- location
      
    } else {
      
      location_code <- percentASFR %>%
        as_tibble() %>% 
        dplyr::filter(.data$name %in% location) %>% 
        select("country_code") %>% 
        unique() %>%
        as.numeric()
      
    }
    
    # find balancing anchor years for interpolation
    # If user requested interpolation for years e.g. 1950 and 1951
    # we will also need a year 1950, i.e. next available year to inerpolate
    # this function finds these anchor years
    anchor_years <- get_bracketing_years(Dates = AsfrDatesIn)
    
    # first prepare tfr data - find country and years
    tfr1 <- tfr %>%
      dplyr::filter(.data$country_code %in% location_code) %>%
      pivot_longer(-c("country_code", "name"),
                   names_to  = "year",
                   values_to = "tfr") %>%
      mutate("year" = suppressWarnings(parse_number(.data$year))) %>% 
      dplyr::filter(.data$year %in% anchor_years)
    
    # now calculate 5-year asfr
    asfr <- percentASFR %>%
      dplyr::filter(.data$country_code %in% location_code) %>%
      pivot_longer(-c("country_code", "name", "age"),
                   names_to  = "year",
                   values_to = "pasfr") %>% 
      mutate("year" = parse_number(.data$year)) %>%
      dplyr::filter(.data$year %in% anchor_years) %>% 
      # join TFR so we can convert percentage-of-TFR to ASFR
      left_join(tfr1, by = c("country_code", "name", "year")) %>%
      # create asfr
      # Convert pasfr (percent of TFR) to absolute ASFR for each single age.
      # Important: we divide pasfr by sum(pasfr)to ensure normalization.
      # this one is for 5 years so no division by 5
      mutate(
        "asfr" = (.data$pasfr / sum(.data$pasfr)) * .data$tfr,
        "age"  = parse_number(.data$age),
        .by  = c("country_code", "name", "year")
      ) %>% 
      select(-c("pasfr", "tfr"))
    
    # now we need female population to calculate births in future
    popage <- popF %>% 
      dplyr::filter(.data$country_code %in% location_code) %>%
      mutate("age" = parse_number(.data$age)) %>% 
      dplyr::filter(.data$age %in% unique(asfr$age)) %>%
      pivot_longer(-c("country_code", "name", "age"),
                   names_to  = "year",
                   values_to = "pop",
                   names_transform = list("year" = as.integer)) %>% 
      dplyr::filter(.data$year %in% anchor_years) %>% 
      # population is per 1000
      mutate("pop" = .data$pop * 1000)
    
    age1 <- range(unique(popage$age))
    
    # now we calculate births and graduate them with pclm using pop as offset
    asfr1 <- asfr %>%
      full_join(popage, by = c("country_code", "name", "year", "age")) %>% 
      mutate("births5" = .data$asfr * .data$pop) %>%
      group_nest(.data$country_code, .data$name, .data$year) %>% 
      mutate(data1 = map(data, ~ graduate_pclm(Value  = .x$births5,
                                               Age    = .x$age,
                                               OAG    = FALSE,
                                               # new upper limit i.e if data is 45 then 49
                                               OAnew  = (max(.x$age) + 4),
                                               offset = .x$pop)%>% 
                           as_tibble() %>% 
                           mutate(age = min(age1):(max(age1) + 4)) %>% 
                           rename("asfr" = "value"))) %>% 
      select(-c("data")) %>%
      unnest("data1")
    
    # now we have options for 5 and single year output
    if(output == "single") { 
      
      age1 <- unique(asfr1$age)
      
      # we need to interpolate between adjacent time
      out <- asfr1 %>% 
        pivot_wider(names_from  = "year",
                    values_from = "asfr") %>% 
        select(-"age") %>% 
        group_nest(.data$country_code, .data$name) %>%
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
        unnest("data") %>% 
        mutate("age" = age1) %>% 
        select(-c("country_code", "name")) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    }
    
    # OUTPUT OPTION: 5-year aggregated ASFR
    # For 5-year groups we use original asfr data and interpolate time
    
    if(output == "5-year") {
      
      age5 <- unique(asfr$age)
      
      out <- asfr %>% 
        # wide format
        pivot_wider(names_from  = "year", 
                    values_from = "asfr") %>%
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
        unnest("data") %>%
        mutate("age" = age5) %>% 
        select(-c("country_code", "name")) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    }
    
  } else { 
    
    # environment
    env <- new.env()
    
    # download mx1dt data from the latest package version
    data("percentASFR1dt", package = latest_wpp, envir = env)
    data("tfr1",           package = latest_wpp, envir = env)
    
    percentASFR1dt <- env$percentASFR1dt
    tfr1           <- env$tfr1
    # ------------------------------------------------------ #
    # find location code from the provided location
    # if location is misspelled return NA
    
    if(is.numeric(location)) {
      
      location_code <- location
      
    } else {
      
      location_code <- percentASFR1dt %>%
        as_tibble() %>% 
        dplyr::filter(.data$name %in% location) %>% 
        select("country_code") %>% 
        unique() %>%
        as.numeric()
      
    }
    
    agen <- unique(percentASFR1dt$age)
    
    # -------------------------------
    # Prepare TFR series filtered by country and year
    # -------------------------------
    tfr <- tfr1 %>% 
      as_tibble() %>% 
      pivot_longer(-c("country_code", "name"),
                   names_to  = "year",
                   values_to = "tfr") %>% 
      dplyr::filter(.data$country_code %in% location_code,
                    .data$year < parse_number(latest_wpp) + 1) %>% 
      mutate("year" = as.integer(.data$year))
    
    # -------------------------------
    # OUTPUT OPTION: single-year ASFR
    # -------------------------------
    # For single-year ASFR we use percentASFR1dt (pasfr is percent by single-year age)
    # and convert percentages to rates using TFR. After constructing wide matrices
    # we interpolate the columns (years) to requested AsfrDatesIn.
    
    if(output == "single") { 
      
      out <- percentASFR1dt %>%
        dplyr::filter(.data$country_code %in% location_code,
                      .data$year <= max(tfr$year)) %>%
        # join TFR so we can convert percentage-of-TFR to ASFR
        left_join(tfr, by = c("country_code", "name", "year")) %>%
        # create asfr
        # Convert pasfr (percent of TFR) to absolute ASFR for each single age.
        # Important: earlier code divided pasfr by sum(pasfr) to ensure normalization.
        mutate("asfr" = (.data$pasfr / sum(.data$pasfr)) * .data$tfr,
               .by = c("country_code", "name", "year")) %>%
        select(-c("pasfr", "tfr")) %>% 
        # wide format
        pivot_wider(names_from  = "year", 
                    values_from = "asfr") %>%
        select(-"age") %>%
        # group rows by country/name in case there are multiple territories
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
        unnest("data") %>%
        mutate("age" = agen) %>% 
        select(-c("country_code", "name")) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    }
    
    # OUTPUT OPTION: 5-year aggregated ASFR
    # For 5-year groups we aggregate percentASFR into 5-year bins first (by summing
    # normalized percentages), then multiply aggregated percentages by TFR to obtain
    # ASFR for each 5-year group, and then interpolate across years.
    if(output == "5-year") {
      
      # age5: the set of 5-year group labels derived from single-year ages
      age5 <- unique((agen %/% 5) * 5)
      
      out <- percentASFR1dt %>% 
        dplyr::filter(.data$country_code %in% location_code,
                      .data$year <= max(tfr$year)) %>%
        # create asfr
        # Convert individual single-year percentages to normalized shares
        mutate("pasfr" = (.data$pasfr / sum(.data$pasfr)),
               "age"   = (.data$age %/% 5) * 5,
               .by = c("country_code", "name", "year")) %>% 
        # sum normalized shares across the 5-year bins
        summarise("pasfr" = sum(.data$pasfr), 
                  .by = c("country_code", "name", "year", "age")) %>%
        # attach TFR and convert to ASFR for 5-year groups
        left_join(tfr, by = c("country_code", "name", "year")) %>% 
        mutate("asfr" = .data$pasfr * .data$tfr) %>%
        select(-c("pasfr", "tfr")) %>% 
        # wide format
        pivot_wider(names_from  = "year", 
                    values_from = "asfr") %>%
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
        unnest("data") %>%
        mutate("age" = age5) %>% 
        select(-c("country_code", "name")) %>% 
        column_to_rownames("age") %>% 
        as.matrix()
      
    } 
  }
  
  # Check and Assign user provided ages if they comply with requirements
  if(prt_msg & length(Age) != nrow(out)) { 
    
    warning("The ages you provided do not match with the requested data. Assigning ages automatically from the wpp data")
    
  } else { 
    
    message("Using the provided Age vector for rownames of resulting nLx matrix")
    
  }
  
  return(out)
  
}



# #' Extract SRB estimates from WPP2019
# #' @description We use the `WPP2019_births` dataset from `DemoToolsData` for the sex ratio at birth. Births from WPP 2019 were graduates to single year totals.
# #' @param SRB sex ratio at birth. Either `NULL`, a scalar to assume constant, or a vector of length 3, assumed.
# #' @param location UN Pop Div `LocName` or `LocID`
# #' @param DatesOut numeric vector of three decimal dates produced by `basepop_ive()`
# #' @param verbose logical, shall we send optional messages to the console?
# #' @return numeric vector with three SRB estimates
# #' @export
# #' @importFrom stats setNames
# 
# 
# downloadSRB <- function(SRB, location, DatesOut, verbose = TRUE){
#   
#   
#   
#   if (!is.null(SRB)) {
#     if (length(SRB) > 3) stop("SRB can only accept three dates at maximum")
#     
#     rep_times <- 3 - length(SRB)
#     SRB <- c(SRB, rep(SRB, times = rep_times))
#     return(stats::setNames(SRB[1:3], DatesOut))
#   }
# 
# 
#   if (length(DatesOut) > 3) stop("SRB can only accept three dates at maximum")
#   WPP2019_births <- DemoToolsData::WPP2019_births
#   SRB_default <- round((1 - .4886) / .4886, 3)
# 
#   if (! is_Loc_available(location)) {
#     if (verbose) {
#       cat(paste(location, "not available in DemoToolsData::WPP2019_births\n"))
#       cat(paste("Assuming SRB to be", SRB_default, "\n"))
#     }
# 
#     return(stats::setNames(rep(SRB_default, 3), DatesOut))
#   }
# 
#   if (verbose){
#     cat(paste0("\nbirths not provided. Downloading births for ", loc_message(location), ", for years between ", round(DatesOut[1], 1), " and ", round(DatesOut[length(DatesOut)], 1), "\n"))
#   }
#   LocID <- get_LocID(location)
#   ind <- WPP2019_births$LocID == LocID &
#        WPP2019_births$Year %in% floor(DatesOut)
#   years_srb <- WPP2019_births[ind, "Year", drop = TRUE]
#   SRB <- stats::setNames(WPP2019_births[ind, "SRB", drop = TRUE], years_srb)
# 
#   if (length(SRB) == 0) return(stats::setNames(rep(SRB_default, 3), DatesOut))
# 
#   DatesOut <- floor(DatesOut)
#   yrs_present <- DatesOut %in% years_srb
#   if (any(!yrs_present)) {
#     yrs_not_present <- mean(SRB[as.character(DatesOut[yrs_present])])
#     yrs_not_present <- stats::setNames(rep(yrs_not_present, sum(!yrs_present)), DatesOut[!yrs_present])
#     SRB <- c(SRB, yrs_not_present)
#   }
# 
#   SRB <- SRB[order(as.numeric(names(SRB)))]
#   SRB
# }

# Author: Tim, Rustam
#' Updated Download or Process Sex Ratio at Birth (SRB)
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
#' @references
#' Keyfitz, N. and Flieger, W. (1971). *Population: Facts and Methods of Demography*. San Francisco: W. H. Freeman.
#'
#' The default proportion female at birth (\code{0.4886}), corresponding to \code{SRB = 1.0486}, is a widely used demographic constant commonly attributed to Keyfitz & Flieger (1971). This explains why the default value is not the rounded \code{1.05}.
#' @importFrom dplyr select rename reframe mutate
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom readr parse_number
#' @importFrom rlang .data
#' @importFrom magrittr %>% 
#' @importFrom stats approx setNames
#' @importFrom utils data
#'
#' @examples
#' # Example 1: Download SRB data for France
#' # lets try SRB 
#' SRB <- downloadSRB(SRB      = NULL,
#'                    location = "Argentina",
#'                    DatesOut = 1950:1952,
#'                    verbose  = TRUE)
#'
#' # another setup for Spain
#' SRB2 <- downloadSRB(SRB      = NULL,
#'                     location = "Spain",
#'                     DatesOut = 2000:2030,
#'                     verbose  = TRUE)
#'
# With user provided SRB.
#'SRB_user <- downloadSRB(SRB = SRB[1],
#'                        DatesOut = 1950)
#'
#' @export

downloadSRB <- function(SRB      = NULL,
                        location = NULL, 
                        DatesOut = NULL, 
                        verbose  = TRUE) {
  
  # 0.4886 = proportion female at birth from Keyfitz & Flieger (1971)
  SRB_default <- round((1 - 0.4886) / 0.4886, 3)
  
  # Check DatesOut
  if(length(DatesOut) < 1 & is.null(SRB)) {
    
    stop("DatesOut must contain at least one date.")
    
  }  
  
  # If SRB provided directly
  # ------------------------------------------------------------------ #
  # CASE 1: User-supplied SRB values (no download required)
  # ------------------------------------------------------------------ #
  if(!is.null(SRB)) {
    
    if(is.null(DatesOut)) {
      
      DatesOut <- names2age(SRB)
      
    }
    
    SRB <- setNames(rep(SRB, length.out = length(DatesOut)), DatesOut)
    
    return(SRB)
    
  }
  
  # ------------------------------------------------------------------ #
  # CASE 2: Retrieve SRB from latest installed WPP package
  # ------------------------------------------------------------------ #
  if(is.null(SRB)) { 
    
    installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
    
    if(length(installed_wpp) == 0) {
      
      stop("No WPP package installed. Please use check_and_load_latest_wpp() function to install the package.")
      
    }
    
    latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
    
    if(parse_number(latest_wpp) < 2022) {
      
      warning("No single-year SRB is available in WPP versions earlier than wpp2022.\n",
              "The 5-year SRB values are interpolated into single years assuming uniform distribution.\n",
              "Please update the wpp package if you want the exact values.")
      
      # Load SRB data from the selected WPP package
      env <- new.env()
      data("sexRatio", package = latest_wpp, envir = env)
      sexRatio1 <- env$sexRatio
      
    } else { 
      
      env <- new.env()
      # Load SRB data from the selected WPP package
      data("sexRatio1", package = latest_wpp, envir = env)
      sexRatio1 <- env$sexRatio1
      
    }
    
    # Determine country code (numeric or character location input)
    if(is.numeric(location)) {
      
      location_code <- location
      
    } else {
      
      location_code <- sexRatio1 %>%
        as_tibble() %>%
        dplyr::filter(.data$name %in% location) %>%
        select("country_code") %>%
        unique() %>%
        as.numeric()
    }
    
    # Handle missing or invalid locations gracefully
    if(is.na(location_code)) {
      
      if(verbose) {  
        
        message(paste0(location, " not available in wpp. Using default SRB: ", SRB_default))
        
      }
      
      return(setNames(rep(SRB_default, length(DatesOut)), DatesOut))
    }
    
    # Notify user of download progress
    if(verbose) {
      
      message(paste0("Downloading SRB for ", 
                     location, 
                     " for years ", 
                     round(min(DatesOut), 1), 
                     " to ", 
                     round(max(DatesOut), 1)
      ))
      
    }
    
    # ------------------------------------------------------------------ #
    # Extract and format SRB data for the specified country and years
    # ------------------------------------------------------------------ #
    dt <- sexRatio1 %>%
      as_tibble() %>%
      dplyr::filter(.data$country_code == location_code) %>%
      pivot_longer(-c("country_code", "name"), 
                   names_to  = "year", 
                   values_to = "SRB") %>%
      mutate("year" = parse_number(.data$year)) %>%
      dplyr::filter(.data$year %in% floor(DatesOut))
    
    # in case we need interpolation we assume uniformity
    if(parse_number(latest_wpp) < 2022) {
      
      dt <- dt %>%
        reframe(
          # single years
          "year_out" = seq(min(DatesOut), max(DatesOut)),
          
          # Original year & SRB vectors for interpolation
          "SRB" = approx(
            x      = .data$year,         # original years
            y      = .data$SRB,          # original SRB values
            xout   = seq(min(.data$year), max(.data$year)),
            method = "constant",         # or "linear"
            f      = 0                   # left-continuous
          )$y,
          .by = c("country_code", "name")
        ) %>%
        rename("year" = "year_out")
      
      years_srb <- dt$year
      SRB       <- setNames(dt$SRB, years_srb)
      
      return(SRB)
      
      # in case data originally comes by single years
    } else {
      
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
  }
  
  return(SRB)
  
}

#' extract births from latest wpp version installed on the machine
#' @param births `NULL` or else a vector of births to simply return
#' @param yrs_births Numeric. vector of years to extract
#' @param location Character or Numeric. UN Pop Dov `LocName` or `LocID`
#' @param sex Character. `"male"`, `"female"`, or `"both"`
#' @param verbose Logical, shall we send optional messages to the console?
#' @return vector of births
#' @export
#' @importFrom utils data
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom dplyr mutate select full_join group_by summarise rename group_nest
#' @importFrom tidyr pivot_longer pivot_wider unnest 
#' @importFrom purrr map
#' @importFrom readr parse_number
#' @examples
#' # Example 1: Download SRB data for France
#' # lets try SRB 
#' births <- fetch_wpp_births(births     = NULL,
#' yrs_births = 1960:1970,
#' location   = "Argentina",
#' sex        = "male",
#' verbose    = TRUE)
#'
# With user provided SRB.
#'SRB_user <- fetch_wpp_births(births = births)
#'
#' @export

fetch_wpp_births <- function(births     = NULL,
                             yrs_births = NULL,
                             location   = NULL,
                             sex        = "female",
                             verbose    = TRUE) {
  
  # if births are provided just return them back
  if(!is.null(births)) {
    
    return(births)
    
  }
  
  # if not provided
  if(is.null(births) | length(births) == 0) {
    
    # check that location, sex, and yrs_births are provided by the user
    if(is.null(location) | length(location) == 0) {
      stop(
        "A minimum of location, sex, and yrs_births should be provided by the user to run the function. location is not provided."
      )
      
    }
    
    if(is.null(sex) | length(sex) == 0) {
      stop(
        "A minimum of location, sex, and yrs_births should be provided by the user to run the function. sex is not provided."
      )
      
    }
    
    if(is.null(yrs_births) | length(yrs_births) == 0) {
      stop(
        "A minimum of location, sex, and yrs_births should be provided by the user to run the function. yrs_births is not provided."
      )
      
    }
    
    # now find the package version
    installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
    
    # stop if none
    if(length(installed_wpp) == 0) {
      stop(
        "No wpp package installed. Please use the function chack_and_load_latest_wpp() to install a wpp package version."
      )
    }
    
    # find the lates one
    latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
    
    # Now in any case to obtain births we need asfr and SRB
    # we already have functions for these so use them here
    # asfr
    asfr <- downloadAsfr(
      Asfrmat     = NULL,
      location    = location,
      AsfrDatesIn = yrs_births,
      verbose     = FALSE,
      output      = "single") %>%
      as.data.frame() %>%
      rownames_to_column("age") %>%
      mutate("age" = as.numeric(.data$age))
    
    
    # special case when user only requested 1 year
    if(length(yrs_births) < 2) { 
      
      names(asfr) <- c("age", yrs_births)
      
    }
    
    asfr <- asfr %>%
      pivot_longer(
        -c("age"),
        names_to  = "year",
        values_to = "asfr",
        names_transform = list(year = as.integer)
      )
    
    # srb
    SRB <- downloadSRB(
      SRB      = NULL,
      DatesOut = yrs_births,
      location = location,
      verbose  = FALSE) %>% 
      as.data.frame() %>%
      rownames_to_column("year") %>%
      mutate("year" = as.numeric(.data$year)) %>% 
      rename("srb" = ".")
    
    # now 
    env <- new.env()
    
    if(any(yrs_births < 1950, yrs_births > (parse_number(latest_wpp) + 1))) {
      message(paste0(
        "Careful, choosing beyond range 1950-",
        parse_number(latest_wpp)
      ))
      
    }
    
    if(parse_number(latest_wpp) < 2022) {
      
      if(verbose) {
        cat(
          paste0(
            "\nbirths not provided. Downloading births for ",
            loc_message(location),
            ", gender: ",
            "`",
            sex,
            "`, years: ",
            paste(yrs_births, collapse = ", "),
            "\n"
          )
        )
      }
      
      # get the bracketing years
      anchor_years <- get_bracketing_years(Dates = yrs_births)
      
      warning(
        "No single ages or years are availabe in wpp versions earlier than wpp2022.
            Please update the wpp package. Currently the pclm-graduated data will be provided."
      )
      
      # popF
      data("popF", package = latest_wpp, envir = env)
      
      if(is.numeric(location)) {
        
        location_code <- location
        
      } else {
        
        location_code <- env$popF %>%
          dplyr::filter(.data$name %in% location) %>%
          select("country_code") %>%
          unique() %>%
          as.numeric()
      }
      
      popF <- env$popF %>%
        dplyr::filter(.data$country_code %in% location_code) %>%
        mutate("age" = parse_number(.data$age)) %>%
        select(-c("country_code", "name")) %>%
        pivot_longer(
          -c("age"),
          names_to  = "year",
          values_to = "pop",
          names_transform = list(year = as.integer)) %>%
        dplyr::filter(.data$year %in% anchor_years) %>%
        mutate("pop" = .data$pop * 1000) %>%
        dplyr::filter(!is.na(.data$pop)) %>% 
        group_nest(.data$year) %>%
        # we graduate ages with pclm
        mutate(
          data = map(
            data,
            ~ graduate_pclm(Value = .x$pop, 
                            Age   = .x$age) %>%
              as.data.frame() %>%
              rownames_to_column("age") %>%
              rename("pop" = ".")
          )
        ) %>%
        unnest("data") %>%
        mutate("age" = as.numeric(.data$age)) %>%
        pivot_wider(names_from  = "year", 
                    values_from = "pop") %>%
        select(-"age") %>%
        as.matrix()
      
      # now we use interp() to interpolate years
      pop <- interp(popF,
                    as.numeric(colnames(popF)),
                    as.numeric(yrs_births),
                    extrap = TRUE) %>%
        as_tibble() %>%
        mutate(age = 0:(nrow(popF) - 1))
      
      # special case when user only requested 1 year
      if(length(yrs_births) < 2) { 
        
        names(pop) <- c(yrs_births, "age")
        
      }
      
      pop <- pop %>%
        pivot_longer(
          -c("age"),
          names_to  = "year",
          values_to = "pop",
          names_transform = list("year" = as.integer)
        ) %>%
        dplyr::filter(.data$age %in% unique(asfr$age))
      
    } else {
      
      if(verbose) {
        cat(
          paste0(
            "\nbirths not provided. Downloading births for ",
            loc_message(location),
            ", gender: ",
            "`",
            sex,
            "`, years: ",
            paste(yrs_births, collapse = ", "),
            "\n"
          )
        )
      }
      
      # in latest versions we have single years data
      # so no need for interpolation and graduation
      data("popF1", package = latest_wpp, envir = env)
      
      
      if(is.numeric(location)) {
        
        location_code <- location
        
      } else {
        
        location_code <- env$popF1 %>%
          dplyr::filter(.data$name %in% location) %>%
          select("country_code") %>%
          unique() %>%
          as.numeric()
      }
      
      pop <- env$popF1 %>%
        dplyr::filter(.data$country_code %in% location_code) %>%
        as_tibble() %>%
        select(-c("country_code", "name")) %>%
        pivot_longer(
          -c("age"),
          names_to  = "year",
          values_to = "pop",
          names_transform = list(year = as.integer)) %>%
        dplyr::filter(.data$year %in% unique(asfr$year)) %>%
        mutate("pop" = .data$pop * 1000) %>%
        dplyr::filter(.data$age %in% unique(asfr$age))
      
    }
    
    # now we can calculate births
    births <- pop %>%
      full_join(asfr, by = c("age", "year")) %>%
      mutate("births" = .data$asfr * .data$pop) %>%
      group_by(.data$year) %>%
      summarise(births = sum(.data$births)) %>%
      full_join(SRB, by = c("year"))
    
    if(sex == "both") {
      
      bt        <- births$births
      names(bt) <- births$year
      
    }
    
    if(sex == "male") {
      
      bt        <- births$births
      bt        <- bt * births$srb / (1 + births$srb)
      names(bt) <- births$year
      
    }
    
    if(sex == "female") {
      
      bt        <- births$births
      bt        <- bt / (births$srb + 1)
      names(bt) <- births$year
      
    }
  }
  
  return(bt)
  
}

# fetch_wpp_births <- function(births, yrs_births, location, sex, verbose) {
#   
#   # fetch WPP births if not provided by user
#   if (is.null(births) | length(births) == 0) {
#     
#     installed_wpp <- grep("^wpp\\d{4}$", rownames(installed.packages()), value = TRUE)
#     
#     if(length(installed_wpp) == 0) {
#       
#       stop("No WPP package installed. Please use check_and_load_latest_wpp() function to install the package.")
#       
#     }
#     
#     latest_wpp <- sort(installed_wpp, decreasing = TRUE)[1]
#     
#     env <- new.env()
#     data("sexRatio", package = latest_wpp, envir = env)
#     sexRatio1 <- env$sexRatio
#     
#     
#     # load WPP births
#     #requireNamespace("DemoToolsData", quietly = TRUE)
#     WPP2019_births <- DemoToolsData::WPP2019_births
#     
#     # filter out location and years
#     ind       <- WPP2019_births$LocID == get_LocID(location) & 
#       WPP2019_births$Year %in% yrs_births
#     b_filt    <- WPP2019_births[ind, ]
#     bt        <- b_filt$TBirths
#     SRB       <- b_filt$SRB
#     
#     # extract births depending on sex
#     if (sex == "both")  births  <- bt
#     if (sex == "male")   births  <- bt * SRB / ( 1 + SRB)
#     if (sex == "female") births  <- bt / (SRB + 1)
#     
#     if (verbose){
#       cat(paste0("\nbirths not provided. Downloading births for ", loc_message(location), ", gender: ", "`", sex, "`, years: ",paste(yrs_births,collapse = ", "), "\n"))
#     } 
#     
#     return(births)
#   } else {
#     # if births are provided simply return them
#     return(births)
#     
#     }
#   
#   
# }

# Author: Rustam
#' Check, Install, and Load the Latest Available `wpp` package version
#'
#' This function checks the \href{https://github.com/PPgp}{PPgp GitHub organization}
#' for the most recent year-based `wpp` package (e.g., `wpp2024`, `wpp2022`),
#' compares it with any installed version on the local machine, optionally installs
#' an update, and loads the selected package into the current R session.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Retrieves all `wpp` repositories from the PPgp GitHub account.
#'   \item Identifies the most recent available WPP dataset.
#'   \item Determines whether a `wpp` package is already installed locally.
#'   \item Compares installed and available versions.
#'   \item Prompts the user to update when running interactively, or automatically updates in non-interactive sessions.
#'   \item Loads the chosen WPP package (latest available or existing installation).
#' }
#'
#' @return
#' Invisibly returns the name of the loaded WPP package (e.g., `"wpp2024"`),
#' and loads that package into the current R session.
#'
#' @note
#' This function requires internet access and the packages `jsonlite`, `httr`,
#' and `remotes`. The user may be prompted to confirm installation during
#' interactive sessions.
#'
#' @examples
#' # Requires internet and may install/update packages
#' check_and_load_latest_wpp()
#'
#' @importFrom jsonlite fromJSON
#' @importFrom httr GET content
#' @importFrom utils installed.packages compareVersion packageVersion
#' @importFrom remotes install_github
#' @export

check_and_load_latest_wpp <- function() {
  
  # ---------------------------------------------------------------------------- #
  # Step 1: Get all wpp-like repositories from GitHub
  repos     <- jsonlite::fromJSON("https://api.github.com/users/PPgp/repos?per_page=100")
  wpp_repos <- repos$name[grepl("^wpp\\d{4}$", repos$name)]
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 2: Extract years and determine the latest
  wpp_years <- as.integer(sub("wpp", "", wpp_repos))
  
  if(length(wpp_years) == 0) {
    
    stop(
      "No WPP packages were found in the PPgp GitHub organization.\n",
      "Expected repositories named like 'wpp2024', 'wpp2022', etc."
    )    
    
  }
  
  latest_year <- max(wpp_years)
  latest_pkg  <- paste0("wpp", latest_year)
  
  message("Latest available wpp package: ", latest_pkg)
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 3: Check installed wpp versions
  installed_pkgs <- utils::installed.packages()[, "Package"]
  installed_wpp  <- installed_pkgs[grepl("^wpp\\d{4}$", installed_pkgs)]
  
  if(length(installed_wpp) == 0) {
    
    message("No wpp package currently installed.")
    installed_version <- NA_character_
    
  } else {
    
    installed_wpp_years <- as.integer(sub("wpp", "", installed_wpp))
    current_pkg         <- installed_wpp[which.max(installed_wpp_years)]
    installed_version   <- as.character(packageVersion(current_pkg))
    
    message("Installed wpp package: ", 
            current_pkg, 
            " (v", installed_version, ")")
    
  }
  
  # ---------------------------------------------------------------------------- #
  # Step 4: Get DESCRIPTION for latest version
  desc_url <- paste0("https://raw.githubusercontent.com/PPgp/", 
                     latest_pkg, 
                     "/main/DESCRIPTION")
  
  desc <- tryCatch(
    httr::content(httr::GET(desc_url), "text", encoding = "UTF-8"),
    error = function(e) stop("Could not read DESCRIPTION for ", latest_pkg)
  )
  
  latest_version <- sub("^Version:\\s+", 
                        "", 
                        grep("^Version:", strsplit(desc, "\n")[[1]], value = TRUE))
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 5: Compare versions
  needs_update <- is.na(installed_version) ||
    utils::compareVersion(installed_version, latest_version) < 0
  
  if(needs_update) {
    
    message("Newer version available: ", latest_pkg, " (v", latest_version, ")")
    
    # Is session interactive?
    if(interactive()) {
      
      ans <- readline(prompt = "Do you want to update to the latest version? [y/n]: ")
      ans <- tolower(ans)
      
    } else {
      
      # If session is non-interactive, assume answer "yes".
      message("Non-interactive session detected -> proceeding with update.")
      ans <- "y"
      
    }
    # -------------------------------
    
    if(ans == "y") {
      
      message("Installing ", latest_pkg, "...")
      remotes::install_github(paste0("PPgp/", latest_pkg), upgrade = "always")
      message("Updated to ", latest_pkg, " (v", latest_version, ")")
      
    } else {
      
      message("Update skipped. Using currently installed version.")
      
    }
    
  } else {
    
    message("Installed version is already up-to-date.")
    ans <- "n"   # ensure defined
    
  }
  # ---------------------------------------------------------------------------- #
  
  # ---------------------------------------------------------------------------- #
  # Step 6: Load the appropriate version
  final_pkg <- ifelse(needs_update && ans == "y", latest_pkg, current_pkg)
  
    library(final_pkg, character.only = TRUE, quietly = TRUE)
  
  return(final_pkg)
}


interp_coh_download_mortality <- function(location, sex, date1, date2, OAnew = 100, verbose){
  
  . <- NULL
  
  date1      <- dec.date(date1)
  date2      <- dec.date(date2)
  
  year1      <- floor(date1) + 1
  year2      <- floor(date2)
  
  year_seq   <- year1:year2
  
  dates_out  <- c(dec.date(date1), year_seq)
  if (verbose){
    cat(paste0("\nlxMat not provided. Downloading lxMat for ", loc_message(location), ", gender: ", "`", sex, "`, for years between ", round(date1, 1), " and ", round(date2, 1), "\n"))
  } 
  
  PX <- suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
                                locations = location,
                                sex = sex)) %>%
    lapply(function(X){
      X[,c("year","x","mx")]
    }) %>%
    lapply(lt_a2s_chunk, OAnew = OAnew) %>%
    lapply(function(X){
      #1 - X$nqx
      lt_id_Ll_S(X$nLx, X$lx, X$Age, X$AgeInt, N = 1)
    }) %>%
    do.call("cbind",.)
  
  
  dimnames(PX)   <- list(0:OAnew, dates_out)
  
  PX[PX > 1]     <- 1
  # discount first and last periods.
  
  f1             <- diff(dates_out)[1]
  f2             <- date2 - floor(date2)
  
  # assume linear px change within age class
  PX[, 1]        <- PX[, 1] ^f1
  PX[,ncol(PX)]  <- PX[, ncol(PX)] ^f2
  
  PX
}

loc_message <- function(location){
  cds     <- DemoToolsData::WPP_codes
  if (fertestr::is_LocID(location)){
    LocName   <- get_LocName(location)
    LocID     <- location
  } else {
    LocID     <- get_LocID(location)
    LocName   <- location
  }
  paste0(LocName," (LocID = ",LocID,")")
  
}

get_LocID <- function(location){
  if (fertestr::is_LocID(location)){
    return(location)
  } else {
    cds     <- DemoToolsData::WPP_codes
    ind       <- cds$LocName == location
    if (!any(ind)){
      stop("requested LocName not found")
    }
    LocID     <- cds[ind,"LocID"] %>% c()
    return(LocID)
  }
}

get_LocName <- function(location){
  if (fertestr::is_LocID(location)){
    cds         <- DemoToolsData::WPP_codes
    ind         <- cds$LocID == location
    if (!any(ind)){
      stop("requested LocID not found")
    }
    LocName     <- cds[ind,"LocName"] %>% c()
    return(LocName)
  } else {
    return(location)
  }
}

is_Loc_available <- function(location){
  isID   <- fertestr::is_LocID(location)
  cds <- DemoToolsData::WPP_codes
  if (isID){
    out <- location %in% cds$LocID
  } else {
    out <- location %in% cds$LocName
  }
  out
}

# function I use for interpolation ancor dates assignment
# in case user requested it to be done within some 5-year anchor points 
# like 1950, 1955. Imagine user requested 1950, 1951
# then interpolation shold consider 2 points 1950, 1955
# similarly, request is 1950, 1956, 1969, then points are 1950, 1955, 1960, 1965, 1970
get_bracketing_years <- function(Dates = NULL) {

  # set up range of available years
  years5 <- seq(from = 1950, to = 2100, by = 5)
  
  # range of input years
  ymin <- min(Dates)
  ymax <- max(Dates)
  
  # find all 5-year anchor years covering this range
  anchors <- years5[years5 >= range(Dates)[1] & years5 <= range(Dates)[2]]
  
  # extend by one anchor below and one above if they exist
  lower_anchor <- max(years5[years5 <= ymin])
  upper_anchor <- min(years5[years5 >= ymax])
  
  # If there is only ONE anchor (e.g., Dates = 2000), 
  # add the NEXT available anchor year.
  if(lower_anchor == upper_anchor) {
    
    upper_anchor <- years5[which(years5 == lower_anchor) + 1]
    
  }
  
  # return unique sorted years
  return(sort(unique(c(lower_anchor, anchors, upper_anchor))))
  
}