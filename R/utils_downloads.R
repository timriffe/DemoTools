
# These utils might be used by basepop, interp_coh, OPAG, mig_resid*,
# and potentially others.

#' Extract Lx estimates from WPP2019
#' @description We use the `FetchLifeTableWpp2019` function of the `fertestr` to extract `Lx` from `wpp2019`, interpolated to an exact date.
#' @param nLx either `NULL` or a numeric vector of lifetable exposure. If it's the second then we just pass it back.
#' @param location UN Pop Div `LocName` or `LocID`
#' @param gender `"male"`, `"female"`, or `"both"`
#' @param nLxDatesIn numeric vector of three decimal dates produced by (or passed through) `basepop_five()`
#'
#' @return numeric matrix of `nLx` with `length(nLxDatesIn)` and abrdiged ages in rows.
#' @export
#'
downloadnLx <- function(nLx, location, gender, nLxDatesIn) {
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
    
    if (is.null(location)) stop("You need to provide a location to download the data for nLx")
    
    if (verbose) {
      cat(paste0("Downloading nLx data for ", location, ", years ", paste(nLxDatesIn,collapse=", "), ", gender ", gender), sep = "\n")
    }
    
    . <- NULL
    
    ind_invalidyear <- which(nLxDatesIn < 1955)
    if (length(ind_invalidyear) != 0) {
      invalid_yrs <- paste0(nLxDatesIn[ind_invalidyear], collapse = ", ")
      cat("nLxDate(s)", invalid_yrs, "is/are below 1955. Capping at 1955\n")
      nLxDatesIn[ind_invalidyear] <- 1955
    }
    
    nLx <-
      lapply(nLxDatesIn, function(x) {
        fertestr::FetchLifeTableWpp2019(location, x, gender)$Lx
      }) %>%
      do.call("cbind", .) %>%
      as.matrix()
    
    colnames(nLx) <- nLxDatesIn
    n             <- nrow(nLx)
    Age           <- c(0,1,seq(5,(n-2)*5,by=5))
    rownames(nLx) <- Age
    return(nLx)
  }
}

downloadAsfr <- function(Asfrmat, location = NULL, AsfrDatesIn) {
  requireNamespace("fertestr", quietly = TRUE)
  verbose <- getOption("basepop_verbose", TRUE)
  
  if (!is.null(Asfrmat)) {
    # TR: can we assume colnames are AsfrDatesIn ?
    return(Asfrmat)
  }
  
  if (is.null(location)) stop("You need to provide a location to download the data for Asfrmat")
  
  ind_invalidyear <- which(AsfrDatesIn < 1955)
  if (length(ind_invalidyear) != 0) {
    invalid_yrs <- paste0(AsfrDatesIn[ind_invalidyear], collapse = ", ")
    cat("AsfrDate(s)", invalid_yrs, "is/are below 1955. Capping at 1955\n")
    AsfrDatesIn[ind_invalidyear] <- 1955
  }
  
  if (verbose) {
    cat(paste0("Downloading Asfr data for ", 
               loc_message(location), 
               ", years ", 
               paste0(AsfrDatesIn),collapse=", "), sep = "\n")
  }
  
  tmp <-
    lapply(AsfrDatesIn, function(x) {
      res        <- fertestr::FetchFertilityWpp2019(location, x)["asfr"]
      names(res) <- NULL
      as.matrix(res)[2:nrow(res), , drop = FALSE]
    })
  
  Asfrmat           <- do.call(cbind, tmp)
  colnames(Asfrmat) <- AsfrDatesIn
  Asfrmat
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
#'
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

