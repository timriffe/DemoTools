
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
  
  tmp <-
    lapply(AsfrDatesIn, function(x) {
      
      if (verbose) {
        cat(paste0("Downloading Asfr data for ", location, ", year ", x), sep = "\n")
      }
      
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
#' @param verbose Whether to print messages. Default set to TRUE.
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
  
  if (!location %in% WPP2019_births$LocName) {
    if (verbose) {
      cat(paste(location, "not available in WPP LocName list\n"))
      cat(paste("Assuming SRB to be", SRB_default, "\n"))
    }
    
    return(stats::setNames(rep(SRB_default, 3), DatesOut))
  }
  
  ind <- WPP2019_births$LocName == location & WPP2019_births$Year %in% floor(DatesOut)
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
