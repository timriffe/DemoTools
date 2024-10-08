#' Net migration by age for an UN family
#' @description Given a total net migration,
#' calculate the net migration age schedule based on the Rogers and Castro formula for UN families.
#' @param NM numeric. Total net migration to distribute between ages and sex.
#' @param family character. Could be "Family", "Female Labor", "Male Labor".
#' @param Single logical. Results by simple age. Default `TRUE`.
#' Typically from pre-working age and working age parts of in Roger-Castro formula.
#' @param OAnew The age from which to group all ages into an open ended age group.
#' By default it is set to 100, so it groups all ages up to 120, which is the
#' maximum age.
#' @export
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @return List with
#' \itemize{
#'   \item{params_RC} {data.frame. Roger-Castro parameters in a data.frame. Same as `mig_un_params` data.}
#'   \item{net_migr} {data.frame. Net migrants by age and sex for the chosen family.}
#' }
#' @examples
#' # 10000 net migrants, comparing two possible families
#' nm1 <- mig_un_fam(NM = 10000, family = "Male Labor", OAnew = 100)
#' nm2 <- mig_un_fam(NM = 10000, family = "Family", OAnew = 100)
#' # See the female profile in for these models:
#' \dontrun{
#' plot(nm1$net_migr$age[nm1$net_migr$sex=="Female"],
#'      nm1$net_migr$nm[nm1$net_migr$sex=="Female"],
#'      xlab="Age",ylab="nm",ylim=c(0,300))
#' points(nm2$net_migr$age[nm2$net_migr$sex=="Female"],
#'        nm2$net_migr$nm[nm2$net_migr$sex=="Female"], col=2)
#' }
mig_un_fam <- function(NM, family, Single = TRUE, OAnew = 100){

  # TR added for global binding warnings
  sex             <- NULL
  age             <- NULL
  .               <- NULL

  mig_un_families <- DemoTools::mig_un_families
  mig_un_params   <- DemoTools::mig_un_params

  mig_sign        <- ifelse(NM < 0, "Emigration", "Inmigration")

  # get asked
  ind         <- mig_un_params$family == family &
                 mig_un_params$mig_sign == mig_sign
  this_params <- mig_un_params[ind,   c("family","sex","param","median")]

  # TR: not priority, but it is also the case that we can do all this with only params
  # see commented-out code below for how to estimate 'family' from params.

  ind         <- mig_un_families$family == family &
                 mig_un_families$mig_sign == mig_sign
  this_family <- mig_un_families[ind,  c("family","sex","age","prop")]

  # get exact 1
  this_family$prop <- this_family$prop + this_family$prop/sum(this_family$prop) * (1-sum(this_family$prop))

  # results
  this_family$nm   <- this_family$prop * NM
  this_family$prop <- NULL

  # Group by family and sex and group ages according to the open
  # age group defined in OAnew.
  this_family <- as.data.table(this_family)
  this_family <- this_family[, .(nm = groupOAG(nm, age, OAnew = OAnew), age = 0:OAnew), by = list(family, sex)]

  # single age
  if(!Single){
    nm              <- NULL
    this_family$age <- trunc(this_family$age/5)*5
    this_family     <- setDT(this_family)[order(sex,age), .(nm=sum(nm)),
                                          by=.(family, age, sex)] %>% as.data.frame()
    }

  # out
  list(net_migr = this_family,
       params_RC = this_params)
}


# data construction -------------------------------------------------------

## library(devtools)
## library(tidyverse)
## load_all()

## # families from UN
## UN_flies <- readxl::read_excel("~/Downloads/UNPD_Migration Age Patterns-Lookup.xlsx",
##                    skip = 3, col_names = T) %>%
##                     rename(Type=1, Age=2) %>%
##                     gather(Sex,Prop,-Age,-Type) %>%
##                     mutate(Prop = Prop/100000)

## # no retirement or old-age pattern
## UN_flies %>% ggplot() +
##               geom_line(aes(Age, Prop, col=Type)) +
##               facet_grid(~Sex) + coord_flip()

## # one case - OK
## db <- UN_flies %>% dplyr::filter(Type == "Female Labor Emigration", Sex == "Male")
## a <- graduate(abs(db$Prop),db$Age,method = "sprague")
## b <- graduate(abs(db$Prop),db$Age,method = "beers(ord)")
## c <- graduate(abs(db$Prop),db$Age,method = "uniform")
## sum(a);sum(b);sum(c);sum(db$Prop)
## plot(db$Age,db$Prop/5,t="s",ylim=c(-.015,.015))
## lines(0:80,-a,col=2)
## lines(0:80,-b,col=3)
## lines(0:80,-c,t="s",col=4)
## res <- mig_estimate_rc(0:120,as.numeric(a),
##                        pre_working_age = TRUE,
##                        working_age = TRUE,
##                        retirement = FALSE,
##                        post_retirement = FALSE)

## lines(0:80, -res[["fit_df"]]$median, col = "violet")
## pars <- res$pars_df$median
## pars <- list(a1 =pars[1],  alpha1 = pars[3],
##              a2 = pars[2], alpha2 = pars[4], mu2 = pars[7], lambda2 = pars[6],
##              c = pars[5])

## ages <- 0:120
## mx_RC <- mig_calculate_rc(ages = ages, pars = pars)
## lines(0:80, -mx_RC, col = "black")
## sum(UN_flies %>% dplyr::filter(Type == "Female Labor Emigration", Sex == "Male") %>% pull(Prop))
## sum(-res[["fit_df"]]$median)

## # fit RC params
## .=NULL
## UN_params <- UN_flies %>% split(list(UN_flies$Sex,UN_flies$Type))
## UN_params <- lapply(names(UN_params),
##                     function(X,M){
##                       x = M[[X]]
##                       x_grad <- data.frame(mx = as.numeric(graduate(abs(x$Prop),x$Age,method = "sprague")),
##                                            Age = 0:max(x$Age))
##                       res <- mig_estimate_rc(x_grad$Age, x_grad$mx,
##                                              pre_working_age = TRUE,
##                                              working_age = TRUE,
##                                              retirement = FALSE,
##                                              post_retirement = FALSE)
##                       params <- res$pars_df
##                       params$Type = unique(x$Type)
##                       params$Sex = unique(x$Sex)
##                       params
##                     }, M = UN_params) %>%
##               do.call("rbind",.)

## # test gof
## UN_estimates <- UN_params %>%
##                   split(list(UN_params$Type,UN_params$Sex))

## UN_estimates <-  lapply(names(UN_estimates),
##                    function(X,M){
##                      x = M[[X]]
##                      pars <- pull(x[,"median"])
##                      params <- c(a1 = pars[1], alpha1 = pars[3],
##                                  a2 = pars[2], alpha2 = pars[4], mu2 = pars[7], lambda2 = pars[6],
##                                   c = pars[5])
##                      ages <- 0:120
##                      out <- data.frame(Type = unique(x$Type),
##                                        Sex = unique(x$Sex),
##                                        Age = ages,
##                                        Prop = mig_calculate_rc(ages, params))
##                      out$Prop <- ifelse(stringr::str_detect(out$Type,"Emigration"),-out$Prop,out$Prop)
##                      out
##                    }, M = UN_estimates) %>%
##                   do.call("rbind",.)

## UN_estimates %>% ggplot() +
##   geom_line(aes(Age, Prop, col=Type)) +
##   facet_grid(~Sex) + coord_flip()

## tolerance_admited <- .005
## test_that("lc w lim data works", {
##   # total
##   expect_equal(
##     UN_flies %>% arrange(Type,Sex) %>%
##       group_by(Type,Sex) %>%
##       summarise(Prop = sum(Prop)) %>% pull(Prop),
##     UN_estimates %>% arrange(Type,Sex) %>%
##       mutate(Age = trunc(Age/5)*5) %>%
##       group_by(Type,Sex) %>%
##       summarise(Prop = sum(Prop)) %>% pull(Prop),
##     tolerance = tolerance_admited)
##   # by age
##   expect_equal(
##     UN_flies %>% arrange(Type,Sex,Age) %>% pull(Prop),
##     UN_estimates %>% arrange(Type,Sex,Age) %>%
##       mutate(Age = trunc(Age/5)*5) %>%
##       group_by(Type,Sex,Age) %>%
##       summarise(Prop = sum(Prop)) %>% pull(Prop),
##     tolerance = tolerance_admited)
##   })

## # save data
## UN_params$family    <- trimws(gsub("Emigration|Immigration", "", UN_params$Type))
## UN_params$mig_sign  <- ifelse(stringr::str_detect(UN_params$Type,"Emigration"),
##                              "Emigration","Inmigration")
## UN_params$param     <- rep(c("a1","a2","alpha1","alpha2","c","lambda2","mu2"),
##                          length(unique(UN_params$Type))*2)
## UN_estimates$family    <- trimws(gsub("Emigration|Immigration", "", UN_estimates$Type))
## UN_estimates$mig_sign  <- ifelse(stringr::str_detect(UN_estimates$Type,"Emigration"),
##                               "Emigration","Inmigration")
## mig_un_params          <- UN_params %>% select(family, sex=Sex, mig_sign, param, median)
## mig_un_families        <- UN_estimates %>% select(family, sex=Sex, mig_sign, age=Age, prop=Prop)
## usethis::use_data(mig_un_params, overwrite = TRUE)
## usethis::use_data(mig_un_families, overwrite = TRUE)
