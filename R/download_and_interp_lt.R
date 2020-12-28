
# This script does nothing yet, still in development, deciding how to 
# graduate abridged lifetables

# This is a temporary script to hold a utility function for
# interp_coh()

# goal will be to fill a mortality surface between two censuses.
# args should be date1, date2, country.

interp_coh_download_mortality <- function(country, sex, date1, date2){
  library(lubridate)
  country <- "France"
  sex <- "male"
  
  date1 <- "1994-07-01"
  date2 <- "2003-05-01"
  
  date1 <- as_date(date1)
  date2 <- as_date(date2)
  
  year1 <- year(date1)
  year2 <- year(date2)
  
  year_seq <- year1:year2
  
  dates_out <- c(dec.date(date1),year_seq[-1],dec.date(date2))
  years_grab <- unique(year_seq - year_seq %% 5)
  
  dx <-lapply(dates_out,fertestr::FetchLifeTableWpp2019,
         locations = country, 
         sex = sex) %>% 
    lapply("[[","dx") %>% 
    dplyr::bind_cols() %>% 
    as.matrix() %>% 
    '['(,1)
  dx %>% 
    graduate_mono(Age=x) %>% 
    lt_id_d_q() %>% 
    plot(x=x, type='l', log = 'y')
  plot(0:100,graduate_mono(dx,Age=x))
  x <- c(0,1,seq(5,100,by=5))
  mx1 <- splinefun(cumsum(mx)~x, method = "monoH.FC")(0:100)
 plot(0:100,mx1)
 plot(exp(diff(mx1)))
  res <- fertestr::FetchLifeTableWpp2019(country, dates_out, sex)$mx
  
}

