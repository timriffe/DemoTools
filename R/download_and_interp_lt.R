
# This script does nothing yet, still in development, deciding how to 
# graduate abridged lifetables

# This is a temporary script to hold a utility function for
# interp_coh()

# goal will be to fill a mortality surface between two censuses.
# args should be date1, date2, country.



lt_a2s_chunk <- function(chunk, OAnew, ...){
  ndx <- chunk$dx
  nLx <- chunk$Lx
  Age <- chunk$x
  
  lt_abridged2single(ndx = ndx, 
                     nLx = nLx, 
                     Age = Age, 
                     OAnew = OAnew, 
                     ...)
}

interp_coh_download_mortality <- function(country, sex, date1, date2, OAnew = 100){

  date1      <- dec.date(date1)
  date2      <- dec.date(date2)
  
  year1      <- floor(date1) + 1
  year2      <- floor(date2)
  
  year_seq   <- year1:year2
  
  dates_out  <- c(dec.date(date1), year_seq)

  
  PX <- suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
                               locations = country, 
                               sex = sex)) %>% 
                          lapply(function(X){
                            X[,c("year","x","dx","Lx")]
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

# lxMat <-suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
#                              locations = country, 
#                              sex = sex) %>% 
#                         lapply("[[","lx") %>% 
#                         dplyr::bind_cols() %>% 
#                         as.matrix())

interp_coh_lxMat_pxt <- function(lxMat, 
                                 dates_lx, 
                                 age_lx, 
                                 date1, 
                                 date2, 
                                 OAnew, ...){
  # TR: this is a temp functin, a stop-gap. Some redundant code with
  # interp_coh_download_mortality(), which it itself temporary.
  # the age graduation will move to lt_abridged2single() as soon as it's
  # fixed.
  date1      <- dec.date(date1)
  date2      <- dec.date(date2)
  
  year1      <- floor(date1) + 1
  year2      <- floor(date2)
  
  year_seq   <- year1:year2
  
  dates_out  <- c(dec.date(date1),year_seq)
  
  # get ndx andnLx from lt_abridged()
  
  a1  <- 0:OAnew
  qx1 <- matrix(ncol = ncol(lxMat),
                nrow = length(a1), 
                dimnames = list(a1,
                                dates_lx))
  for (i in 1:ncol(lxMat)){
    LTA     <- lt_abridged(Age = age_lx, 
                           lx = lxMat[, i], 
                           OAnew = OAnew, 
                           radix = 1e6,
                           ...)
    
    LT1     <- lt_abridged2single(ndx = LTA$ndx, 
                                  nLx = LTA$nLx, 
                                  Age = LTA$Age, 
                                  OAnew = OAnew,
                                  ...)
    qx1[, i] <- LT1$nqx
  }
  
  # We do linear interpolation of the logit-transformed qx.
  logit_qx  <- log(qx1 / (1 - qx1))
  
  logit_qx_interp     <- 
                       interp(
                         popmat = logit_qx, 
                         datesIn = dates_lx,
                         datesOut = dates_out,
                         rule = 2) 
  # transform back
  QX            <- exp(logit_qx_interp) / (1 + exp(logit_qx_interp))
  
  
  QX[nrow(QX)]  <- 1
  
  
  f1            <- diff(dates_out)[1]
  f2            <- date2 - floor(date2)
  
  # assume linear px change within age class
  PX            <- 1 - QX
  PX[,1]        <- PX[, 1] ^f1
  PX[,ncol(PX)] <- PX[, ncol(PX)] ^f2

  
  PX
}
