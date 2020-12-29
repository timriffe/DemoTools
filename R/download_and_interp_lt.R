
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
  mono_fun  <- stats::splinefun(CDF~x)
  
  xout      <- 0:100
  sqxout    <- sqrt(xout)
  CDF1      <- mono_fun(sqxout)
  dx1       <- diff(c(0,CDF1)) 
  qx1       <- lt_id_d_q(dx1)
  
  tibble(Age = xout, Date = rep(Date, length(xout)), qx = qx1)
}


interp_coh_download_mortality <- function(country, sex, date1, date2){

  date1      <- dec.date(date1)
  date2      <- dec.date(date2)
  
  year1      <- floor(date1) + 1
  year2      <- floor(date2)
  
  year_seq   <- year1:year2
  
  dates_out  <- c(dec.date(date1),year_seq)
 
  DX <-lapply(dates_out,fertestr::FetchLifeTableWpp2019,
         locations = country, 
         sex = sex) %>% 
    lapply("[[","dx") %>% 
    dplyr::bind_cols() %>% 
    as.matrix()
  
  dimnames(DX) <- list(c(0,1,seq(5,100,by=5)),
                  dates_out)

  
  DXlong <- reshape2::melt(DX, varnames = c("Age","Date"), value.name = "dx")
  
  QX <-
    DXlong %>% 
    dplyr::group_by(Date) %>% 
    dplyr::do(graduate_dx_qx_chunk(chunk = .data)) %>% 
    dplyr::ungroup() %>% 
    reshape2::acast(Age~Date, value.var = "qx")
  
  # discount first and last periods.
  
  f1    <- diff(dates_out)[1]
  f2    <- date2 - floor(date2)
  
  QX[,1] <- QX[,1] * f1
  QX[,2] <- QX[,2] * f2

  PX <- 1 - QX

  PX
}

# interp_coh_download_mortality("France","male","1971-07-01","1978-07-01")
