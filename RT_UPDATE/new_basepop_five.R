
basepop_five <- function(location     = NULL,
                         refDate      = NULL,
                         Age          = NULL,
                         Females_five = NULL,
                         Males_five   = NULL,
                         nLxFemale    = NULL,
                         nLxMale      = NULL,
                         nLxDatesIn   = NULL,
                         AsfrMat      = NULL,
                         AsfrDatesIn  = NULL,
                         SRB          = NULL,
                         SRBDatesIn   = NULL,
                         radix        = NULL,
                         verbose      = TRUE,
                         ...) {
  
  options(basepop_verbose = verbose)
  on.exit(options(basepop_verbose = NULL))
  refDate <- dec.date(refDate)
  
  if(!is.null(Age)) {
    
    stopifnot(is_abridged(Age))
    stopifnot(length(Age) == length(Females_five))
    
  } else {
    
    if(!is.null(names(Females_five))) {
      
      Age <- names2age(Females_five)
      
    } else {
      
      if(verbose) {
        
        cat("Assuming age groups are in standard abridged intervals")
        
      }
      
      Age <- inferAgeIntAbr(Females_five)
      
    }
  }
  
  if(is.null(nLxDatesIn)) {
    
    nLxDatesIn <- refDate - c(0.5, 7.5)
    
    if(verbose) {
      
      cat(paste0(
        "Assuming the two prior dates for the nLx matrix to be: ",
        paste0(nLxDatesIn, collapse = ", ")
      ),
      sep = "\n")
      
    }
  }
  
  if(is.null(AsfrDatesIn)) {
    
    AsfrDatesIn <- refDate - c(0.5, 7.5)
    
    if(verbose){
      
      cat(paste0(
        "Assuming the two prior dates for the Asfr matrix to be: ",
        paste0(AsfrDatesIn, collapse = ", ")
      ),
      sep = "\n")
      
    }
  }
  names(Females_five) <- Age
  names(Males_five)   <- Age
  
  nLxFemale <- download_nLx(
    nLx        = nLxFemale,
    location   = location,
    gender     = "female",
    nLxDatesIn = sort(nLxDatesIn),
    output     = "abridged",
    radix      = 1
  )
  
  nLxMale <- download_nLx(
    nLx        = nLxFemale,
    location   = location,
    gender     = "male",
    nLxDatesIn = sort(nLxDatesIn),
    output     = "abridged",
    radix      = 1
  )
  
  
  if(is.null(radix)) {
    
    radix <- lt_infer_radix_from_1L0(nLxMale[1, 1])
    
    if(verbose) {
      cat(
        paste0(
          "Setting radix to value of lx: ",
          radix,
          ". Can be overwritten with the `radix` argument"
        ),
        sep = "\n"
      )
    }
  }
  
  # not ok
  AsfrMat <- download_Asfr(Asfrmat     = Asfrmat,
                           location    = location,
                           AsfrDatesIn = sort(AsfrDatesIn),
                           output      = "5-year")
  
  
  DatesOut   <- refDate - c(0.5, 2.5, 7.5)
  SRBDatesIn <- if(!is.null(SRBDatesIn)) {
    
    SRBDatesIn
    
  } else { 
    
    DatesOut
    
  }
  
  SRB <- downloadSRB(SRB       = NULL,
                     location  = location,
                     DatesOut  = sort(SRBDatesIn),
                     verbose  = verbose)
  
  AllArgs <- as.list(environment())
  ArgsCheck(AllArgs)
  lower_bound <- abs(min(nLxDatesIn) - min(DatesOut))
  upper_bound <- abs(max(nLxDatesIn) - max(DatesOut))
  if (lower_bound > 5 || upper_bound > 5) {
    stop(
      "nLxDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates"
    )
  }
  
  nLxf <- interp(nLxFemale, datesIn = nLxDatesIn, datesOut = DatesOut)
  nLxm <- interp(nLxMale,   datesIn = nLxDatesIn, datesOut = DatesOut)
  
  lower_bound <- abs(min(AsfrDatesIn) - min(DatesOut))
  upper_bound <- abs(max(AsfrDatesIn) - max(DatesOut))
  
  if(lower_bound > 5 || upper_bound > 5) {
    stop(
      "AsfrDatesIn implies an extrapolation of > 5 years to achieve the needed reference dates"
    )
  }
  
  Asfr                <- interp(AsfrMat, datesIn = AsfrDatesIn, datesOut = DatesOut)
  Asfr                <- Asfr[-c(1, ncol(Asfr)), ]
  ages_15_55          <- as.character(seq(15, 55, by = 5))
  ages_20_55          <- ages_15_55[-1]
  ages_15_50          <- ages_15_55[-9]
  ages_20_50          <- ages_15_55[-c(1, 9)]
  ages_15_45          <- ages_15_55[-c(8, 9)]
  ages_20_45          <- ages_15_55[-c(1, 8, 9)]
  ages_15_40          <- ages_15_55[-c(7, 8, 9)]
  FMiddleages         <- Females_five[ages_15_55]
  Ft_minus_5          <- FMiddleages[ages_20_55] * nLxf[ages_15_50, 2] / nLxf[ages_20_55, 2]
  names(Ft_minus_5)   <- ages_15_50
  Ft_minus_10         <- Ft_minus_5[ages_20_50] * nLxf[ages_15_45, 3] / nLxf[ages_20_50, 3]
  names(Ft_minus_10)  <- ages_15_45
  Ft_minus_.5         <- FMiddleages[ages_15_45] * 0.9 + Ft_minus_5[ages_15_45] * 0.1
  Ft_minus_2.5        <- FMiddleages[ages_15_45] * 0.5 + Ft_minus_5[ages_15_45] * 0.5
  Ft_minus_7.5        <- Ft_minus_5[ages_15_45] * 0.5 + Ft_minus_10[ages_15_45] * 0.5
  fExpos              <- cbind(Ft_minus_.5, Ft_minus_2.5, Ft_minus_7.5)
  Bt                  <- colSums(fExpos * Asfr)
  Males_five_out      <- Males_five
  Females_five_out    <- Females_five
  PF                  <- 1 / (SRB + 1)
  Females_five_out[1] <- Bt[1] * PF[1] * nLxf[1, 1] / radix
  Males_five_out[1]   <- Bt[1] * (1 - PF[1]) * nLxm[1, 1] / radix
  Females_five_out[2] <- Bt[2] * PF[2] * 5 * sum(nLxf[1:2, 2]) / (radix * 5) - Females_five_out[1]
  Males_five_out[2]   <- Bt[2] * (1 - PF[2]) * 5 * sum(nLxm[1:2, 2]) / (radix * 5) - Males_five_out[1]
  Females_five_out[3] <- Bt[3] * PF[3] * 5 * sum(nLxf[1:2, 3]) / (radix * 5) * nLxf[3, 2] / sum(nLxf[1:2, 2])
  Males_five_out[3] <- Bt[3] * (1 - PF[3]) * 5 * sum(nLxm[1:2, 3]) / (radix * 5) * nLxm[3, 2] / sum(nLxm[1:2, 2])
  
  return(
    list(
      Females_adjusted = Females_five_out,
      Males_adjusted   = Males_five_out,
      Females_five     = Females_five,
      Males_five       = Males_five,
      nLxf             = nLxf,
      nLxm             = nLxm,
      Asfr             = Asfr,
      Exposure_female  = fExpos,
      Bt               = Bt,
      SRB              = SRB,
      Age              = Age
    ))
}