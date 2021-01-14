
# This script does nothing yet, still in development, deciding how to 
# graduate abridged lifetables

# This is a temporary script to hold a utility function for
# interp_coh()

# goal will be to fill a mortality surface between two censuses.
# args should be date1, date2, country.

interp_coh_graduate_dx_qx_chunk <- function(chunk){
  
  # This is a temporary internal method: graduate
  # dx using mono spline on CDF and sqrt(age) (works OK)
  # then convert to qx using identity.
  x         <- sqrt(chunk$Age)
  Date      <- chunk$Date[1]
  CDF       <- cumsum(chunk$dx)
  mono_fun  <- stats::splinefun(CDF~x, method = "monoH.FC")
  
  xout      <- 0:100
  sqxout    <- sqrt(xout)
  CDF1      <- mono_fun(sqxout)
  dx1       <- diff(c(0,CDF1)) 
  qx1       <- lt_id_d_q(dx1)
  
  tibble(Age = xout, Date = rep(Date, length(xout)), qx = qx1)
}


interp_coh_download_mortality <- function(country, sex, date1, date2, OAnew = 100){

  date1      <- dec.date(date1)
  date2      <- dec.date(date2)
  
  year1      <- floor(date1) + 1
  year2      <- floor(date2)
  
  year_seq   <- year1:year2
  
  dates_out  <- c(dec.date(date1),year_seq)
 
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
  
  PX <-suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
                               locations = country, 
                               sex = sex)) %>% 
                          lapply(function(X){
                            X[,c("year","x","dx","Lx")]
                          }) %>% 
                          dplyr::bind_rows() %>% 
    group_by(year) %>% 
    do(lt_a2s_chunk(chunk = .data, OAnew = OAnew)) %>% 
    mutate(px = 1-nqx) %>% 
    reshape2::acast("Age~year", value.var = "px")
  
  PX[PX > 1] <- 1
  # discount first and last periods.
  
  f1    <- diff(dates_out)[1]
  f2    <- date2 - floor(date2)
  
  # assume linear px change within age class
  PX[, 1]        <- PX[,1] ^f1
  PX[,ncol(PX)] <- PX[, ncol(PX)] ^f2

  PX
}

# lxMat <-suppressMessages(lapply(dates_out,fertestr::FetchLifeTableWpp2019,
#                              locations = country, 
#                              sex = sex) %>% 
#                         lapply("[[","lx") %>% 
#                         dplyr::bind_cols() %>% 
#                         as.matrix())

interp_coh_lxMat_pxt <- function(lxMat, dates_lx, Age_lx, date1, date2, ...){
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
  
  lxfull     <- interp(
                  popmat = lxMat, 
                  datesIn = dates_lx,
                  datesOut = dates_out,
                  rule = 2) 
  rownames(lxfull) <- Age_lx
  lxlong <- reshape2::melt(lxfull, 
                           varnames = c("Age","Date"), 
                           value.var = "lx")
  QX <-
    lxlong %>% 
    dplyr::group_by(.data$Date) %>% 
    dplyr::mutate(dx = lt_id_l_d(.data$lx)) %>% 
    dplyr::do(interp_coh_graduate_dx_qx_chunk(chunk = .data)) %>% 
    dplyr::ungroup() %>% 
    reshape2::acast(Age~Date, value.var = "qx")
  QX[QX < 0] <- 0
    
  
  f1    <- diff(dates_out)[1]
  f2    <- date2 - floor(date2)
  
  # assume linear px change within age class
  PX            <- 1 - QX
  PX[,1]        <- PX[,1] ^f1
  PX[,ncol(PX)] <- PX[,ncol(PX)] ^f2
  # QX[, 1] <- QX[, 1] ^ (1/f1)
  # QX[, ncol(QX)] <- QX[,ncol(QX)] ^ (1/f2)
  
  # PX <- 1 - QX
  
  PX
}
