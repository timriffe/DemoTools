#===============================================================================
# 2019-12-03 -- DemoTools
# Wrap the latest HMD data -- only works locally on ik machine
# Ilya Kashnitsky, ilya.kashnitsky@gmail.com
#===============================================================================


library(tidyverse)
library(fs)

# this one is specific for my (Ilya) machine 
# a local copy of HMD downloaded on 2019-11-28
hmdpath <- fs::as_fs_path("~/data/hmd/")

# wrap data.table::fread to read HMD files
fread_hmd <- function(x) x %>% data.table::fread(skip = 2, na.strings = ".")

# generalize to purrr::map_df the chose directory
fread_hmd_dir <- function(thedir) {
    thedir %>% 
        dir_ls() %>% 
        map_df(fread_hmd, .id = "country") %>% 
        janitor::clean_names() %>% 
        mutate(
            country = country %>% 
                str_remove(thedir %>% path_real()) %>% 
                str_remove("\\..*") %>% 
                str_remove("\\/")
        ) %>% 
        rename(Lx = lx_2, Tx = tx)
}



# read in all HMD 5x5 life tables  ----------------------------------------

lt_both_5x5 <- path(hmdpath, "lt_both", "bltper_5x5") %>% fread_hmd_dir()
lt_female_5x5 <- path(hmdpath, "lt_female", "fltper_5x5") %>% fread_hmd_dir()
lt_male_5x5 <- path(hmdpath, "lt_male", "mltper_5x5") %>% fread_hmd_dir()

hmd_lt_5x5 <-
    bind_rows(lt_both_5x5, lt_female_5x5, lt_male_5x5, .id = "sex") %>%
    mutate(
        sex = sex %>% as_factor %>% 
            lvls_revalue(c("b", "f", "m")) %>% 
            paste
    ) %>% 
    drop_na()


# train wilmoth models, sex specific --------------------------------------

library(MortalityEstimate)

# Fit Log-quadratic model
x <- c(0,1, seq(5, 110, by = 5))
fitted_logquad_b <- wilmoth(x = x, LT = hmd_lt_5x5 %>% filter(sex == "b"))
fitted_logquad_f <- wilmoth(x = x, LT = hmd_lt_5x5 %>% filter(sex == "f"))
fitted_logquad_m <- wilmoth(x = x, LT = hmd_lt_5x5 %>% filter(sex == "m"))


# shrink the size of the fitted model object
fitted_logquad_b$input$LT <- NULL
fitted_logquad_f$input$LT <- NULL
fitted_logquad_m$input$LT <- NULL

fitted_logquad_b$residuals <- NULL
fitted_logquad_f$residuals <- NULL
fitted_logquad_m$residuals <- NULL



# save the fitted models to be used in the package -----------------------

usethis::use_data(fitted_logquad_b, overwrite = T, compress = "xz")
usethis::use_data(fitted_logquad_f, overwrite = T, compress = "gzip")
usethis::use_data(fitted_logquad_m, overwrite = T, compress = "bzip2")


# benchmark different compression read-in speed ---------------------------

microbenchmark::microbenchmark(
    xz = load("data/fitted_logquad_b.rda"),
    gzip = load("data/fitted_logquad_f.rda"),
    bzip2 = load("data/fitted_logquad_m.rda")
)

# Unit: milliseconds
# expr       min        lq      mean    median        uq      max neval
# xz 11.070890 12.761848 13.591440 13.279367 13.723383 18.35390   100
# gzip  1.512547  1.662455  1.914923  1.739654  1.840855 13.81812   100
# bzip2 20.602232 20.876711 21.568214 21.055588 21.499222 25.56252   100



# save all with "gzip" compression ----------------------------------------

usethis::use_data(fitted_logquad_b, overwrite = T, compress = "gzip")
usethis::use_data(fitted_logquad_f, overwrite = T, compress = "gzip")
usethis::use_data(fitted_logquad_m, overwrite = T, compress = "gzip")

