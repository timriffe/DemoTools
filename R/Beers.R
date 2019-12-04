

###############################################################################
# Based on BeersMSplit.R, whose header was:
########################################################################
## Interpolation of of event data (five year periods)
## fifth/single years by Beers Modified six-term formula
## (Siegel and Swanson, 2004, p. 729)
## This formula applies some smoothing to the interpolant, which is
## recommended for time series of events with some dynamic
## R implementation by Thomas Buettner (21 Oct. 2015)
########################################################################


#' Create the Beers ordinary or modified coefficient matrix
#'
#' @description The resulting coefficient matrix is based on the number of rows in \code{Value}
#' which must be in 5-year age groups (not abridged). The final row may be an open
#' or closed age group, as indicated by the \code{OAG} argument.
#'
#' @inheritParams graduate
#' @param method character. Valid values are \code{"mod"} or \code{"ord"}. Default \code{"mod"}.
#'
#' @details The \code{Value} vector is a placeholder in this case. This function is
#' a utility called by the Beers family of functions, where it is most convenient to just pass
#' in the same matrix being used in those calculations to determine the layout of the coefficient matrix.
#'
#' @references
#' \insertRef{beers1945modified}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @export
#' @examples
#' coefsOA     <- graduate_beers_expand(pop5_mat, OAG = TRUE, method = "mod")
#' coefsclosed <- graduate_beers_expand(pop5_mat, OAG = FALSE, method = "mod")
#' dim(graduate_beers_expand(pop5_mat, TRUE))
#' dim(graduate_beers_expand(pop5_mat, FALSE))
#' coefso      <- graduate_beers_expand(pop5_mat, OAG = TRUE, method = "ord")
#'
#' # how to use (under the hood in beers()
#'
graduate_beers_expand <- function(Value,
                        OAG = FALSE,
                        method = "Mod") {
  method <- tolower(method)
  stopifnot(method %in% c("mod", "ord"))
  
  # nr 5-year age groups
  m      <- length(Value)
  # nr rows in coef mat.
  n      <- m * 5 - ifelse(OAG, 4, 0)
  # number of middle blocks
  MP     <- m - ifelse(OAG, 5, 4)
  
  if (method == "mod") {
    ## Beers Modified Split
    g1g2 <- matrix(
      c( 0.3332, -0.1938, 0.0702, -0.0118, 
         0.0022 , 0.2569, -0.0753, 0.0205, 
         -0.0027, 0.0006 , 0.1903, 0.0216, 
         -0.0146, 0.0032, -0.0005 , 0.1334, 
         0.0969, -0.0351, 0.0059, -0.0011 , 
         0.0862, 0.1506, -0.0410, 0.0054,
         -0.0012 , 0.0486, 0.1831, -0.0329, 
         0.0021, -0.0009 , 0.0203, 0.1955, 
         -0.0123, -0.0031, -0.0004 , 0.0008, 
         0.1893, 0.0193, -0.0097, 0.0003 ,
         -0.0108, 0.1677, 0.0577, -0.0153, 
         0.0007 ,-0.0159, 0.1354, 0.0972, 
         -0.0170, 0.0003
      ),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )
    
    g3 <- matrix(
      c( -0.0160, 0.0973, 0.1321, -0.0121, 
         -0.0013 ,-0.0129, 0.0590, 0.1564, 
         0.0018, -0.0043 ,-0.0085, 0.0260, 
         0.1650, 0.0260, -0.0085 ,-0.0043, 
         0.0018, 0.1564, 0.0590, -0.0129 ,
         -0.0013, -0.0121, 0.1321, 0.0973, 
         -0.0160
      ),
      5,
      5,
      byrow = TRUE
    )
    
    g4g5 <- matrix(
      c( 0.0003, -0.0170, 0.0972, 0.1354, 
         -0.0159 , 0.0007, -0.0153, 0.0577, 
         0.1677, -0.0108 , 0.0003, -0.0097, 
         0.0193, 0.1893, 0.0008 ,-0.0004, 
         -0.0031, -0.0123, 0.1955, 0.0203 ,
         -0.0009, 0.0021, -0.0329, 0.1831, 
         0.0486 ,-0.0012, 0.0054, -0.0410, 
         0.1506, 0.0862 ,-0.0011, 0.0059, 
         -0.0351, 0.0969, 0.1334 ,-0.0005, 
         0.0032, -0.0146, 0.0216, 0.1903 , 
         0.0006, -0.0027, 0.0205, -0.0753, 
         0.2569 , 0.0022, -0.0118, 0.0702, 
         -0.1938, 0.3332
      ),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )
  }
  if (method == "ord") {
    ## Beers Oscillatory (Johnson spreadsheet)
    g1g2 <- matrix(
      c( 0.3333, -0.1636, -0.021, 0.0796, 
         -0.0283, 0.2595, -0.078, 0.013, 
         0.01, -0.0045, 0.1924, 0.0064, 
         0.0184, -0.0256, 0.0084, 0.1329, 
         0.0844, 0.0054, -0.0356, 0.0129, 
         0.0819, 0.1508, -0.0158, -0.0284, 
         0.0115, 0.0404, 0.2, -0.0344, 
         -0.0128, 0.0068, 0.0093, 0.2268, 
         -0.0402, 0.0028, 0.0013,-0.0108, 
         0.2272, -0.0248, 0.0112, -0.0028,
         -0.0198, 0.1992, 0.0172, 0.0072, 
         -0.0038,-0.0191, 0.1468, 0.0822, 
         -0.0084, -0.0015),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )
    
    g3 <- matrix(
      c( -0.0117, 0.0804, 0.157, -0.0284, 
         0.0027,-0.002, 0.016, 0.22, 
         -0.04, 0.006, 0.005, -0.028, 
         0.246, -0.028, 0.005, 0.006, 
         -0.04, 0.22, 0.016, -0.002, 
         0.0027, -0.0284, 0.157, 0.0804, 
         -0.0117),
      5,
      5,
      byrow = TRUE
    )
    
    g4g5 <- matrix(
      c( -0.0015, -0.0084, 0.0822, 0.1468, 
         -0.0191,-0.0038, 0.0072, 0.0172, 
         0.1992, -0.0198,-0.0028, 0.0112, 
         -0.0248, 0.2272, -0.0108, 0.0013, 
         0.0028, -0.0402, 0.2268, 0.0093, 
         0.0068, -0.0128, -0.0344, 0.2, 
         0.0404, 0.0115, -0.0284, -0.0158, 
         0.1508, 0.0819, 0.0129, -0.0356, 
         0.0054, 0.0844, 0.1329, 0.0084, 
         -0.0256, 0.0184, 0.0064, 0.1924,
         -0.0045, 0.01, 0.013, -0.078, 
         0.2595,-0.0283, 0.0796, -0.021, 
         -0.1636, 0.3333
      ),
      nrow = 10,
      ncol = 5,
      byrow = TRUE
    )
  }
  
  ## create a Beers coefficient matrix for 5-year age groups
  bm               <- matrix(0, nrow = n, ncol =  m)
  ## insert upper left block
  bm[1:10, 1:5]    <- g1g2
  
  # determine positions of middle blocks
  rowpos           <-
    matrix(11:((MP * 5) + 10), ncol = 5, byrow = TRUE)
  colpos           <- row(rowpos) + col(rowpos) - 1
  for (i in (1:MP)) {
    # calculate the slices and add middle panels accordingly
    bm[rowpos[i,], colpos[i,]] <- g3
  }
  ## standard format for Beers coefficients
  
  ## insert last two panels
  
  fr                <- nrow(bm) - ifelse(OAG, 10, 9)
  lr                <- fr + 9
  fc                <- ncol(bm) - ifelse(OAG, 5, 4)
  lc                <- fc + 4
  bm[fr:lr, fc:lc]   <- g4g5
  
  
  if (OAG) {
    # preserve open ended age group
    bm[nrow(bm), ncol(bm)]    <- 1
  }
  
  bm
}

#' The ordinary modified Beers splitting methods
#'
#' @description This method offers both ordinary and modified Beers splitting, with an optional \href{https://www.census.gov/data/software/dapps.html}{Demographic Analysis & Population Projection System Software} adjustment \code{johnson} for ages under 10.
#'
#' @inheritParams graduate
#' @param method character. Valid values are \code{"mod"} or \code{"ord"}. Default \code{"mod"}.
#' @param johnson  logical. Whether or not to adjust young ages according to the \href{https://www.census.gov/data/software/dapps.html}{Demographic Analysis & Population Projection System Software} method. Default \code{FALSE}.
#' @details Ages should refer to lower age bounds. \code{Value} must be labelled with ages unless \code{Age} is given separately. There must be at least six 5-year age groups (including the open group, 5 otherwise). If you want the \code{johnson} adjustment then \code{Value} must contain a single-year estimate of the population count in age 0. That means \code{Value} must come either as standard abridged or single age data.
#'
#' If the highest age does not end in a 0 or 5, and \code{OAG == TRUE}, then the open age will be grouped down to the next highest age ending in 0 or 5. If the highest age does not end in a 0 or 5, and \code{OAG == FALSE}, then results extend to single ages covering the entire 5-year age group.
#'
#' @return A numeric vector of single age data.
#' @references
#' \insertRef{beers1945modified}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' \insertRef{stover2008spectrum}{DemoTools}
#' @export
#' @examples
#' p5 <- pop5_mat
#' a5 <- as.integer(rownames(p5))
#' head(p5) # this is the entire matrix
#' p1 <- graduate_beers(p5[,1], Age = a5, OAG = FALSE)
#' head(p1)
#' # note some negatives in high ages
#' tail(p1)
#' sum(p1) - sum(p5[,1])
#'
#' # another case, starting with single ages
#' # note beers() groups ages.
#' Value        <- pop1m_ind
#' Age          <- 0:100
#' names(Value) <- Age
#' ord1 <-  graduate_beers(Value, Age, OAG = TRUE, method = "ord")
#' mod1 <- graduate_beers(Value, Age, OAG = TRUE, method = "mod")
#' \dontrun{
#' plot(Age,Value,
#' ylab = 'The counts', xlab = 'Age groups')
#' lines(Age, ord1, t='l', col='blue')
#' lines(Age, mod1, t = 'l', col ='red', lty =2)
#' legend(80,15000000,
#'       legend = c('Ordinary',
#'                  'Modified'),
#'       col=c('blue', 'red'),
#'       lty = c(1,2))
#' }
#'
#' # notice this negative value. Yuck!
#' tail(mod1)
#' # this replaces ages 90+, guaranteed no negatives.
#' graduate_mono_closeout(Value, Age = Age, pops = mod1, OAG = TRUE)
#' # Note: there are no kludges built into beers() to handle such cases.
#' # these ought to be handled by wrappers as appropriate.
#'
#' # This replicates Johnson_2016_BEERSP.XLS, males
#' M <- c(184499,752124-184499,582662,463534,369976,286946,235867,
#' 		199561,172133,151194,131502,113439,95614,
#' 		78777,60157,40960,21318,25451)
#' Age <- c(0,1,seq(5,80,by=5))
#' Age0        <- 184499
#' johnson     <- graduate_beers(
#' 		         Value = M,
#' 		         Age = Age,
#' 		         OAG = TRUE,
#' 			       method = "ord",
#' 			       johnson = TRUE)
graduate_beers <- function(Value,
                  Age,
                  AgeInt,
                  OAG = TRUE,
                  method = "mod",
                  johnson = FALSE) {
  if (as.character(match.call()[[1]]) == "beers") {
    warning("please use graduate_beers() instead of beers().", call. = FALSE)
  }
  
  if (missing(AgeInt)){
    AgeInt <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  
  punif1       <- graduate_uniform(
                    Value = Value, 
                    AgeInt = AgeInt,
                    Age = Age, 
                    OAG = OAG)

  # this is innocuous if ages are already grouped
  a1           <- as.integer(names(punif1))
  pop5         <- groupAges(
                    punif1, 
                    Age = a1,
                    N = 5,
                    shiftdown = 0)
  # depending on OAG, highest age may shift down.
  a5           <- as.integer(names(pop5))
  # punif1       <- graduate_uniform(
  #                   pop5,
  #                   Age = a5,
  #                   OAG = OAG)
  
  # generate coefficient matrix
  bm                <- graduate_beers_expand(pop5, OAG = OAG, method = method)
  
  # redistribute
  pop1              <- bm %*% pop5
  dim(pop1)         <- NULL
  # can only do the Johnson adjust if ages are single or abridged.
  # cuz we need a separate age 0
  if (johnson & ((min(Age) == 0 & 1 %in% Age))) {
    Age0 <- Value[1]
    pop1 <- johnsonAdjust(Age0 = Age0,
                          pop5 = pop5,
                          pop1 = pop1)
  }
  
  # TR will this fail if Age starts with 1?
  zero              <- min(Age)
  ages              <- zero:(nrow(bm) - 1 + zero)
  names(pop1)       <- ages
  pop1
}

#' @export
#' @rdname graduate_beers
beers <- graduate_beers

#' Adjust ages under 10 using a modification of Beers
#' @description Assuming we have an external estimate of age 0, this method
#' refits to the ordinary Beers single age results, remaining constrained to the
#' original 5-year age groups and smoothly blending into ages greater than 10.
#' @param Age0 numeric. An estimate of age 0.
#' @param pop5 numeric. Matrix of age-period population counts in 5-year age groups.
#' @param pop1 numeric. Matrix of age-period population using Beers ordinary (or some other) method.
#'
#' @return A matrix of single age population estimates.
#' @details This has not been tested using \code{pop1} as generated from other methods, such as
#' the Beers modified, Sprague, or Grabill methods. Called internally by \code{beers()}.
#' @export
#' @references
#' \insertRef{stover2008spectrum}{DemoTools}
johnsonAdjust <- function(Age0, pop5, pop1) {
  # coefficient matrix
  DAPPSmod <- matrix(
    c(
      0.2333, 0.3445, -0.1222, 0.3778,
      -0.1556,-0.0875, 0.3847, -0.0236,
      0.0389, -0.0111,-0.1458,
      0.2708, 0.1458, -0.4167, 0.1667,-0.0833,
      0.1151, 0.2817, -0.6111, 0.2222,
      0.0000, -0.0080, 0.3254, -0.3889,
      0.1111, 0.0458, -0.0613, 0.2637,
      0.1833, -0.1000, 0.0375, -0.0458,
      0.1292, 0.8167, -0.2333
    ),
    ncol = 5,
    byrow = TRUE
  )
  # the composite pop thing
  py          <- c(pop1[2],
                       pop5[1] - (Age0 + pop1[2]),
                       sum(pop1[6:9]),
                       pop1[10:11])
  # gives new ages 2-8
  pynew       <- DAPPSmod %*% py
  dim(pynew)  <- NULL
  # recompose output (still constrained)
  pop1[1]   <- Age0   # keep Age0 est
  pop1[3:9] <- pynew
  pop1
}
