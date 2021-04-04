
# These utils might be used by basepop, interp_coh, OPAG, mig_resid*,
# and potentially others.

#' Extract Lx estimates from WPP2019
#' @description We extract `Lx` from `wpp2019`, interpolated to exact dates. Different methods availables. 
#' A vector of countries can handle, but with an unique sex. Row names are not indicative of countries.
#' @param nLx numeric. either `NULL` or a numeric vector of lifetable exposure. If it's the second then we just pass it back.
#' @param location vector. UN Pop Div `LocName` or `LocID`
#' @param gender character. `"male"`, `"female"`, or `"both"`
#' @param nLxDatesIn numeric. Vector of three decimal dates produced by (or passed through) `basepop_five()`
#' @param method character. Could be `"linear"`, `"exponential"`, or `"power"`
#'
#' @return numeric matrix of `nLx` with `length(nLxDatesIn)` and abrdiged ages in rows.
#' @export
#' @examples
#' # life expectancy calculated from Lx downloaded from WPP19. Using names or codes.
#' Lxs_name <- downloadnLx(nLx=NULL, location = "Argentina",
#'                         gender = "both", nLxDatesIn = 1950:2030)
#' Lxs_code <- downloadnLx(nLx=NULL, location = "32",
#'                         gender = "both", nLxDatesIn = 1950:2030)
#' \dontrun{
#' plot(1950:2030, as.numeric(colSums(Lxs_name)), xlab = "Year", ylab="e0")
#' lines(1950:2030, as.numeric(colSums(Lxs_code)))
#' }
#' # life expectancy for different countries
#' Lxs_countries <- downloadnLx(nLx=NULL, location = c("Argentina","Brazil","Uruguay"),
#' gender = "both", nLxDatesIn = 1950:2025)
#' \dontrun{
#' plot(1950:2025, as.numeric(colSums(Lxs_countries[1:22,])), 
#'      t="l", xlab = "Year", ylab="e0", ylim = c(40,80))
#' lines(1950:2025, as.numeric(colSums(Lxs_countries[23:44,])), col=2)
#' lines(1950:2025, as.numeric(colSums(Lxs_countries[45:64,])), col=3)
#' legend("bottomright",c("Argentina","Brazil","Uruguay"),lty=1,col=1:3)
#' }
downloadnLx <- function(nLx, location, gender, nLxDatesIn, method="linear") {
  
  requireNamespace("fertestr", quietly = TRUE)
  requireNamespace("magrittr", quietly = TRUE)
  verbose <- getOption("basepop_verbose", TRUE)
  
  if (!is.null(nLx)) {
    # TR: ensure colnames passed
    nLx <- as.matrix(nLx)
    colnames(nLx) <- nLxDatesIn
    n             <- nrow(nLx)
    Age           <- c(0,1,seq(5,(n-2)*5,by=5))
    rownames(nLx) <- Age
    return(nLx)
  }
  
  if (is.null(nLx)){
    
    # stop/warnings
    if (is.null(location)){
      stop("You need to provide a location to download the data for nLx")
    } 
    if (!any(is_LocID(location))) {
      location_code <- get_location_code(location)
    }else {
      location_code <- as.integer(location)
    }
      
    if (verbose) {
      cat(paste0("Downloading nLx data for ", location, ", years ", paste(nLxDatesIn,collapse=", "), ", gender ", gender), sep = "\n")
    }
    if(any(nLxDatesIn<1950,nLxDatesIn>2025)){
      cat("Careful, extrapolating beyond range 1950-2025")
    }
    
    # handle sex
    sex_code <- ifelse(tolower(gender) == "both", "b", 
                       ifelse(tolower(gender) == "female", "f", 
                              ifelse(tolower(gender) == "male", "m", NA)))
    Sex_mortlaws <- ifelse(sex_code == "b", "total", tolower(gender))
    stopifnot(`Invalid sex name, please set it to 'both', 'male' or 'female'` = !is.na(sex_code))
    
    # initial data
    lt_wpp19 <-DemoToolsData::WPP2019_lt
    
    # filter and matrix shape
    lt_ctry <- lt_wpp19[lt_wpp19$LocID %in% location_code &
                          lt_wpp19$Sex %in% sex_code,] %>% as.data.frame() %>% 
                reshape(data = ., 
                        direction = "wide", idvar = c("LocID","AgeStart","Sex"), 
                        timevar = "Year", v.names = "mx", drop = c("AgeSpan","lx"))

    # intert/extrap rates and built life tables for each combination location/Sex/Year 
    .<-NULL
    out <- cbind(lt_ctry[,c(1:3)],
                 interp(lt_ctry[,-c(1:3)], 
                        seq(1953,2023,5), as.numeric(nLxDatesIn), 
                        extrap = TRUE, method = method) %>% 
                        as.data.frame() %>% 
                   setNames(as.character(nLxDatesIn))
                 ) %>% 
          split(., list(lt_ctry$LocID, lt_ctry$Sex)) %>% 
          lapply(function(X){
            Age <- X[["AgeStart"]]
            apply(X[,-c(1:3)] %>% 
                    as.data.frame()%>% setNames(as.character(nLxDatesIn)), 2,
                    function(S){
                        MortalityLaws::LifeTable(x = Age,
                                                 mx = S,
                                                 lx0 = 1,
                                                 sex = Sex_mortlaws)$lt$Lx
                  })
          }) %>% 
          do.call("rbind", .)
    
    # combination as rowname
    rownames(out) <- lt_ctry$AgeStart
    
    return(out)
  }
}

#' Extract ASFR estimates from WPP2019
#' @description We extract `ASFRx` from `wpp2019`, interpolated to exact dates. Different methods availables.
#' A vector of countries can handle, but with an unique sex. Row names are not indicative of countries.
#' @param Asfrmat numeric.
#' @param location vector. UN Pop Div `LocName` or `LocID`
#' @param AsfrDatesIn numeric. Vector of decimal dates.
#' @param method character. Could be `"linear"`, `"exponential"`, or `"power"`
#'
#' @return numeric matrix interpolated asfr
#' @export
#' @examples
#' # Total fertility ratio calculated from ASFRx downloaded from WPP19. 
#' # See `downloadnLx` for analogous examples on multiple countries or using codes instead of names. 
#' ASFR_Arg <- downloadAsfr(Asfrmat = NULL, location = "Argentina", AsfrDatesIn = 1950:2025)
#' \dontrun{
#' plot(1950:2025, as.numeric(colSums(ASFR_Arg))*5, xlab = "Year", ylab="TFR", ylim=c(1.5,4), t="l")
#' }
downloadAsfr <- function(Asfrmat, location = NULL, AsfrDatesIn, method="linear") {
  requireNamespace("fertestr", quietly = TRUE)
  verbose <- getOption("basepop_verbose", TRUE)
  
  if (!is.null(Asfrmat)) {
    # TR: can we assume colnames are AsfrDatesIn ?
    return(Asfrmat)
  }
  
  # stop/warnings
  if (is.null(location)){
    stop("You need to provide a location to download the data for Asfrmat")
  } 
  if (!any(is_LocID(location))) {
    location_code <- get_location_code(location)
  }else {
    location_code <- as.integer(location)
  }
  if (verbose) {
    cat(paste0("Downloading ASFR data for ", location, ", years ", paste(AsfrDatesIn,collapse=", ")), sep = "\n")
  }
  if(any(AsfrDatesIn<1950,AsfrDatesIn>2025)){
    cat("Careful, extrapolating beyond range 1950-2025")
  }
  
  # initial data
  asfr_wpp19    <-DemoToolsData::WPP2019_asfr
  
  # spread format
  asfr_ctry     <- asfr_wpp19[asfr_wpp19$LocID %in% location_code,] %>% 
                      as.data.frame() %>% 
                      reshape(direction = "wide", idvar = c("LocID","AgeStart"), 
                              timevar = "Year", v.names = "ASFR")
  
  # interp/extrap
  out <- interp(asfr_ctry[,-c(1:3)], seq(1953,2023,5),
                as.numeric(AsfrDatesIn), 
                extrap = TRUE, method = method) %>% 
                as.data.frame() %>% 
                setNames(as.character(AsfrDatesIn)) %>% 
                as.matrix()
  
  # combination as rowname
  rownames(out) <- asfr_ctry$AgeStart
  
  return(out)
}

#' Extract SRB estimates from WPP2019
#' @description We use the `WPP2019_births` dataset from `DemoToolsData` for the sex ratio at birth. Births from WPP 2019 were graduates to single year totals.
#' @param SRB sex ratio at birth. Either `NULL`, a scalar to assume constant, or a vector of length 3, assumed.
#' @param location UN Pop Div `LocName` or `LocID`
#' @param DatesOut numeric vector of three decimal dates produced by `basepop_ive()`
#' @param verbose logical, shall we send optional messages to the console?
#'
#' @return numeric vector with three SRB estimates
#' @export


downloadSRB <- function(SRB, location, DatesOut, verbose = TRUE){
  
  
  
  if (!is.null(SRB)) {
    if (length(SRB) > 3) stop("SRB can only accept three dates at maximum")
    
    rep_times <- 3 - length(SRB)
    SRB <- c(SRB, rep(SRB, times = rep_times))
    return(stats::setNames(SRB[1:3], DatesOut))
  }
  

  if (length(DatesOut) > 3) stop("SRB can only accept three dates at maximum")
  WPP2019_births <- DemoToolsData::WPP2019_births
  SRB_default <- round((1 - .4886) / .4886, 3)
  
  if (! is_Loc_available(location)) {
    if (verbose) {
      cat(paste(location, "not available in DemoToolsData::WPP2019_births\n"))
      cat(paste("Assuming SRB to be", SRB_default, "\n"))
    }
    
    return(stats::setNames(rep(SRB_default, 3), DatesOut))
  }
  
  if (verbose){
    cat(paste0("\nbirths not provided. Downloading births for ", loc_message(location), ", for years between ", round(DatesOut[1], 1), " and ", round(DatesOut[length(DatesOut)], 1), "\n"))
  } 
  LocID <- get_LocID(location)
  ind <- WPP2019_births$LocID == LocID & 
       WPP2019_births$Year %in% floor(DatesOut)
  years_srb <- WPP2019_births[ind, "Year", drop = TRUE]
  SRB <- stats::setNames(WPP2019_births[ind, "SRB", drop = TRUE], years_srb)
  
  if (length(SRB) == 0) return(stats::setNames(rep(SRB_default, 3), DatesOut))
  
  DatesOut <- floor(DatesOut)
  yrs_present <- DatesOut %in% years_srb
  if (any(!yrs_present)) {
    yrs_not_present <- mean(SRB[as.character(DatesOut[yrs_present])])
    yrs_not_present <- stats::setNames(rep(yrs_not_present, sum(!yrs_present)), DatesOut[!yrs_present])
    SRB <- c(SRB, yrs_not_present)
  }
  
  SRB <- SRB[order(as.numeric(names(SRB)))]
  SRB
}


#' extract births from wpp2019
#' @param births `NULL` or else a vector of births to simply return
#' @param yrs_births vector of years to extract
#' @param location UN Pop Dov `LocName` or `LocID`
#' @param sex `"male"`, `"female"`, or `"both"`
#' @param verbose logical, shall we send optional messages to the console?
#' @return vector of births
#' @export
#' @importFrom fertestr is_LocID
#' @importFrom fertestr get_location_code
fetch_wpp_births <- function(births, yrs_births, location, sex, verbose) {
  
  # fetch WPP births if not provided by user
  if (is.null(births)) {
    
    # load WPP births
    requireNamespace("DemoToolsData", quietly = TRUE)
    WPP2019_births <- DemoToolsData::WPP2019_births
    
    
   
    
    # filter out location and years
    ind       <- WPP2019_births$LocID == get_LocID(location) & 
      WPP2019_births$Year %in% yrs_births
    b_filt    <- WPP2019_births[ind, ]
    bt        <- b_filt$TBirths
    SRB       <- b_filt$SRB
    
    # extract births depending on sex
    if (sex == "both")  births  <- bt
    if (sex == "male")   births  <- bt * SRB / ( 1 + SRB)
    if (sex == "female") births  <- bt / (SRB + 1)
    
    if (verbose){
      cat(paste0("\nbirths not provided. Downloading births for ", loc_message(location), ", gender: ", "`", sex, "`, years: ",paste(yrs_births,collapse = ", "), "\n"))
    } 
  }
  
  births
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
      1 - X$nqx
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
  if (is_LocID(location)){
    LocName   <- get_LocName(location)
    LocID     <- location
  } else {
    LocID     <- get_LocID(location)
    LocName   <- location
  }
  paste0(LocName," (LocID = ",LocID,")")
  
}

get_LocID <- function(location){
  if (is_LocID(location)){
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
  if (is_LocID(location)){
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
  isID   <- is_LocID(location)
  cds <- DemoToolsData::WPP_codes
  if (isID){
    out <- location %in% cds$LocID
  } else {
    out <- location %in% cds$LocName
  }
  out
}

