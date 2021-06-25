# -------------------------------------------------------------------------- #
# Uniform                                                                    #
# -------------------------------------------------------------------------- #

#' Convert arbitrary age groupings into single years of age.
#'
#' @description Uniformly splits aggregate counts in age groups into single year age groups.
#'
#' @inheritParams graduate
#' @param OAvalue Desired width of open age group. See details.
#'
#' @return Numeric vector of counts for single year age groups.
#'
#' @details Assumes that the population is uniformly distributed across each age interval, and that initial age intervals are integers greater than or equal to 1. If \code{AgeInt} is given, its final value is used as the interval for the final age group. If \code{AgeInt} is missing, then \code{Age} must be given, and the open age group is by default preserved \code{OAvalue} rather than split. To instead split the final age group into, e.g., a 5-year age class, either give \code{AgeInt}, *or* give \code{Age}, \code{OAG = TRUE}, and \code{OAvalue = 5}.  `Age` be any age range, it does not need to start at 0.
#'
#' @export
#' @examples
#' MalePop <- c(9544406,7471790,11590109,11881844,11872503,12968350,
#' 		11993151,10033918,14312222,8111523,15311047,6861510,13305117,7454575,
#' 		9015381,10325432,9055588,5519173)
#' Ages <- seq(0, 85, by = 5)
#' graduate_uniform(MalePop, Age = Ages)
graduate_uniform <-
  function(Value,
           Age,
           AgeInt,
           OAG = TRUE,
           OAvalue = 1) {

    if (missing(Age) & missing(AgeInt)) {
      Age      <- names2age(Value)
    }
    if (missing(AgeInt)) {
      # give 1 to final interval to preserve
      AgeInt   <- age2int(Age, OAG = OAG, OAvalue = OAvalue)
    }
    if (missing(Age)) {
      Age      <- cumsum(AgeInt) - AgeInt
    }
    # discount for single
    out        <- rep(Value / AgeInt, times = AgeInt)
    names(out) <- min(Age):(min(Age) + length(out) - 1)
    out
  }

# -------------------------------------------------------------------------- #
# Sprague                                                                    #
# -------------------------------------------------------------------------- #

#' The basic Sprague age-splitting method.
#'
#' @description This method is used to interpolate counts based on the Sprague formula. It is based on the first stage of the Sprague R script prepared by Thomas Buettner and Patrick Gerland, itself based on the description in Siegel and Swanson, 2004, p. 727.
#'
#' @inheritParams graduate
#' @details Ages should refer to lower age bounds, ending in the open age group in the last row (not a closed terminal age). Dimension labeling is necessary. There must be at least six age groups (including the open group). One year of data will work as well, as long as it's given as or coercible to a single-column matrix. This method may produce negative values, most likely in the youngest or oldest ages. This case is dealt with in the \code{graduate()} wrapper function but not in this function.
#'
#' If the highest age does not end in a 0 or 5, and \code{OAG == TRUE}, then the open age will be grouped down to the next highest age ending in 0 or 5. If the highest age does not end in a 0 or 5, and \code{OAG == FALSE}, then results extend to single ages covering the entire 5-year age group.
#'
#' @return Numeric vector of counts split into single ages.
#'
#' @references
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @export
#'
#' @examples
#' head(pop5_mat) # this is the entire matrix
#' a5 <- as.integer(rownames(pop5_mat))
#' # the last value is an open age group, preserve as such:
#' p1   <- graduate_sprague(Value = pop5_mat[,1], Age = a5, OAG = TRUE)
#' head(p1); tail(p1)
#' sum(p1) - sum(pop5_mat[,1])
#'
#' # another case, starting with single ages
#' Age   <- 0:100
#' # notice how this particular case produces a negative value in the last age
#' # before OAG:
#' pops  <- graduate_sprague(Value = pop1m_ind, Age = Age, OAG = TRUE)
#' # the graduate() wrapper deals with this automatically.
#' \dontrun{
#'   plot(seq(0,100,by=5), pop5_mat[,1]/5, type = 's')
#'   lines(0:100,
#'     p1,
#'     lty = 1,
#'     col = "red")
#'
#' }

graduate_sprague <- function(Value,
                             Age,
                             AgeInt,
                             OAG = TRUE) {

  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
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
    shiftdown = 0
  )
  # depending on OAG, highest age may shift down.
  a5           <- as.integer(names(pop5))

  # generate coefficient matrix
  scm          <- graduate_sprague_expand(
    Value = pop5,
    Age = a5,
    OAG = OAG)

  # redistribute
  pop1         <- scm %*% pop5

  dim(pop1)    <- NULL
  # label and return
  names(pop1)  <- min(Age):(length(pop1) - 1)

  # no sense adding closeout behavior here, when it isn't offered
  # in grabill or beers. Better make wrapper with this sugar.

  #	# default closeout with monoCloseout().
  #	# set to FALSE to turn off, write "mono"
  #	if (is.logical(closeout)){
  #		if (!closeout){
  #			return(pop1)
  #		}
  #		closeout                    <- "mono"
  #	}
  #	if (closeout == "mono"){
  #
  #		# note, if OAG = FALSE and Age %% 5 != 0,
  #		# then we need to group popmat to next lowest
  #		# age divisible by 5...
  #		if (nrow(popmat) > nrow(pop1)){
  #			n                          <- nrow(pop1)
  #			popmat[n, ]                <- colSums(popmat[n:nrow(popmat), ,drop = FALSE])
  #			popmat                     <- popmat[1:n, , drop = FALSE]
  #		}
  #
  #		pop1                        <- monoCloseout(
  #				popmat = popmat,
  #				pops = pop1,
  #				OAG = OAG,
  #				pivotAge = pivotAge)
  #	}
  pop1
}

#' Create the Sprague coefficient matrix.
#'
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat} where is assumed that each row of data is a 5-year age group. The final row may be an open or closed age group, as indicated by the \code{OAG} argument.
#'
#' @inheritParams graduate
#'
#' @details The \code{popmat} matrix is really just a placeholder in this case. This function is a utility called by the Sprague family of functions, where it is most convenient to just pass in the same matrix being used in those calculations to determine the layout of the coefficient matrix.
#'
#' @export
#'
#' @references
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' @examples
#' a5           <- as.integer(rownames(pop5_mat))
#' coefsOA      <- graduate_sprague_expand(pop5_mat[,1], Age = a5, OAG = TRUE)
#' coefsclosed  <- graduate_sprague_expand(pop5_mat[,1], Age = a5, OAG = FALSE)
#' dim(coefsOA)
#' dim(coefsclosed)

graduate_sprague_expand    <- function(
  Value,
  Age,
  OAG = TRUE) {
  popmat                         <- as.matrix(Value)

  # figure out ages and years
  Age5                           <- Age
  Age1                           <- min(Age5):max(Age5)

  # nr 5-year age groups
  m                              <- length(Value)
  # nr rows in coef mat.
  n                              <- m * 5 - ifelse(OAG, 4, 0)
  # number of middle blocks
  MP                             <- m - ifelse(OAG, 5, 4)

  # get the split coefficients
  # block for ages 0-9
  
  # TR: 5-5-2021, this assumes ages start at 0...
  g1g2                           <- matrix(
    c( 0.3616, -0.2768, 0.1488, -0.0336,
       0.0000, 0.2640, -0.0960, 0.0400,
       -0.0080, 0.0000, 0.1840, 0.0400,
       -0.0320, 0.0080, 0.0000, 0.1200,
       0.1360, -0.0720, 0.0160, 0.0000,
       0.0704, 0.1968, -0.0848, 0.0176,
       0.0000, 0.0336, 0.2272, -0.0752,
       0.0144, 0.0000, 0.0080, 0.2320,
       -0.0480, 0.0080, 0.0000,-0.0080,
       0.2160, -0.0080, 0.0000, 0.0000,
       -0.0160, 0.1840, 0.0400, -0.0080,
       0.0000,-0.0176, 0.1408, 0.0912,
       -0.0144, 0.0000
    ),
    nrow = 10,
    ncol = 5,
    byrow = TRUE
  )
  # block for middle ages


  g3                             <- matrix(
    c( -0.0128, 0.0848, 0.1504, -0.0240,
       0.0016,-0.0016, 0.0144, 0.2224,
       -0.0416, 0.0064, 0.0064, -0.0336,
       0.2544, -0.0336, 0.0064, 0.0064,
       -0.0416, 0.2224, 0.0144, -0.0016,
       0.0016, -0.0240, 0.1504, 0.0848,
       -0.0128
    ),
    5,
    5,
    byrow = TRUE
  )

  # block prior to closeout
  g4g5                           <- matrix(
    c(0.0000,-0.0144,0.0912,0.1408,
      -0.0176,0.0000,-0.0080,0.0400,
      0.1840,-0.0160,0.0000,0.0000,
      -0.0080,0.2160,-0.0080,0.0000,
      0.0080,-0.0480,0.2320,0.0080,
      0.0000,0.0144,-0.0752,0.2272,
      0.0336,0.0000,0.0176,-0.0848,
      0.1968,0.0704,0.0000,0.0160,
      -0.0720,0.1360,0.1200,0.0000,
      0.0080,-0.0320,0.0400,0.1840,
      0.0000,-0.0080,0.0400,-0.0960,
      0.2640,0.0000,-0.0336,0.1488,
      -0.2768,0.3616
    ),
    nrow = 10,
    ncol = 5,
    byrow = TRUE
  )

  ## create a Sprague coefficient matrix for 5-year age groups
  bm                             <- matrix(0, nrow = n, ncol =  m)
  ## insert upper left block
  bm[1:10, 1:5]                  <- g1g2

  # determine positions of middle blocks
  rowpos                         <-
    matrix(11:((MP * 5) + 10), ncol = 5, byrow = TRUE)
  colpos                         <- row(rowpos) + col(rowpos) - 1
  for (i in (1:MP)) {
    # calculate the slices and add middle panels accordingly
    bm[rowpos[i,], colpos[i,]]   <- g3
  }

  ## insert last two panels

  fr                             <- nrow(bm) - ifelse(OAG, 10, 9)
  lr                             <- fr + 9
  fc                             <- ncol(bm) - ifelse(OAG, 5, 4)
  lc                             <- fc + 4
  bm[fr:lr, fc:lc]               <- g4g5

  if (OAG) {
    # preserve open ended age group
    bm[nrow(bm), ncol(bm)]       <- 1
  }

  bm
}


# -------------------------------------------------------------------------- #
# Grabill                                                                    #
# -------------------------------------------------------------------------- #

#' Create the Grabill coefficient matrix.
#'
#' @description The resulting coefficient matrix is based on the number of rows in \code{popmat} where we assume that each row of data is a 5-year age group and the final row is an open age group to be preserved as such.
#'
#' @inheritParams graduate
#'
#' @details The \code{Value} vector is really just a placeholder in this case. This function is a utility called by the Grabill family of functions, where it is most convenient to just pass in the same matrix being used in those calculations to determine the layout of the coefficient matrix. Note that these coefficients do not constrain population counts to their year totals. This function is called by \code{grabill()}, which ensures matching marginals by 1) blending boundary ages into the Sprague estimated population, and 2) a second constraint on the middle age groups to enforce matching sums.
#'
#' @references
#' \insertRef{shryock1973methods}{DemoTools}
#'
#' @export
#' @examples
#' a5 <- as.integer(rownames(pop5_mat))
#' graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = TRUE)
#' graduate_grabill_expand(pop5_mat[,1], Age = a5, OAG = FALSE)
graduate_grabill_expand <- function(Value, Age, OAG = TRUE) {
  # figure out ages and years
  Age5   <- Age
  Age1   <- min(Age5):max(Age5)

  # nr 5-year age groups
  m      <- length(Value)
  # nr rows in coef mat.
  n      <- m * 5 - ifelse(OAG, 4, 0)
  # number of middle blocks
  MP     <- m - ifelse(OAG, 5, 4)

  # primary grabill coef block
  g3g    <- matrix(
    c( 0.0111, 0.0816, 0.0826, 0.0256,
       -0.0009, 0.0049, 0.0673, 0.0903,
       0.0377, -0.0002, 0.0015, 0.0519,
       0.0932, 0.0519, 0.0015,-0.0002,
       0.0377, 0.0903, 0.0673, 0.0049,
       -0.0009, 0.0256, 0.0826, 0.0816,
       0.0111
    ),
    5,
    5,
    byrow = TRUE
  )
  ## create a Grabill coefficient matrix for 5-year age groups
  gm                <- matrix(0, nrow = n, ncol =  m)

  fr                <- nrow(gm) - ifelse(OAG, 10, 9)
  lr                <- fr + 9
  fc                <- ncol(gm) - ifelse(OAG, 5, 4)
  lc                <- fc + 4

  # ----------------------------------------------------------
  # Note: for the boundary ages we keep shuffling in g3g, the same grabill
  # coefs. The columns on the boundaries will NOT sum to 1. These coefs are
  # used just for the firs pass, then results blend into the Sprague boundary
  # estimates.
  # ----------------------------------------------------------
  # the young age coefficients
  g1g2g              <- matrix(0, nrow = 10, ncol = 5)
  g1g2g[1:5, 1:3]    <- g3g[, 3:5]
  g1g2g[6:10, 1:4]   <- g3g[, 2:5]
  # the old age coefficients
  g4g5g              <- matrix(0, nrow = 10, ncol = 5)
  g4g5g[1:5, 2:5]    <- g3g[, 1:4]
  g4g5g[6:10, 3:5]   <- g3g[, 1:3]

  gm[1:10, 1:5]      <- g1g2g
  gm[fr:lr, fc:lc]   <- g4g5g


  # determine positions of middle blocks
  rowpos             <- matrix(11:((MP * 5) + 10), ncol = 5, byrow = TRUE)
  colpos             <- row(rowpos) + col(rowpos) - 1
  for (i in (1:MP)) {
    # calculate the slices and add middle panels accordingly
    gm[rowpos[i, ], colpos[i,]] <- g3g
  }

  if (OAG) {
    # preserve open ended age group
    gm[nrow(gm), ncol(gm)]    <- 1
  }

  # return coefficient matrix
  gm
}

#' The basic Grabill age-splitting method
#'
#' @description This method uses Grabill's redistribution of middle ages and blends into
#' Sprague estimated single-age population counts for the first and final ten ages. Open age groups are preserved, as are annual totals.
#'
#' @inheritParams graduate
#' @details  Dimension labeling is necessary. There must be at least six age groups (including the open group). One year of data will work as well, as long as it's given as a single-column matrix. Data may be given in either single or grouped ages. If the highest age does not end in a 0 or 5, and \code{OAG == TRUE}, then the open age will be grouped down to the next highest age ending in 0 or 5. If the highest age does not end in a 0 or 5, and \code{OAG == FALSE}, then results extend to single ages covering the entire 5-year age group.
#'
#' @return numeric vector in single ages.
#'
#' @references
#' \insertRef{shryock1973methods}{DemoTools}
#'
#' @export
#'
#' @examples
#' a5 <- as.integer(rownames(pop5_mat))
#' p5 <- pop5_mat[,1]
#' p1g <- graduate_grabill(Value = p5, Age = a5, OAG = TRUE)
#' sum(p1g) - sum(p5)
#' p1s <- graduate_sprague(p5, Age = a5, OAG = TRUE)
#' \dontrun{
#' plot(seq(0,100,by=5),p5[,1]/5,type = "s", col = "gray", xlab = "Age", ylab = "Count")
#' lines(0:100, p1g, col = "red", lwd = 2)
#' lines(0:100, p1s, col = "blue", lty = 2, lwd =2)
#' legend("topright",
#'		lty = c(1,1,2),
#'		col = c("gray","red","blue"),
#'		lwd = c(1,2,1),
#'		legend = c("grouped","Grabill", "Sprague"))
#' }
#'
#' # also works for single ages:
#' grab1 <- graduate_grabill(Value = pop1m_ind, Age = 0:100)
#' \dontrun{
#' plot(0:100, pop1m_ind)
#' lines(0:100, grab1)
#' }
graduate_grabill <- function(
  Value,
  Age,
  AgeInt,
  OAG = TRUE) {

  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
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
    shiftdown = 0
  )
  # depending on OAG, highest age may shift down.
  a5           <- as.integer(names(pop5))
  # punif1       <- graduate_uniform(
  #                   pop5,
  #                   Age = a5,
  #                   OAG = OAG)


  # get coefficient matrices for Sprague and Grabill
  scmg              <- graduate_grabill_expand(pop5, Age = a5, OAG = OAG)
  scm               <- graduate_sprague_expand(pop5, Age = a5, OAG = OAG)

  # split pop counts
  pops              <- scm %*% pop5
  popg              <- scmg %*% pop5

  # ---------------------------------------------
  # now we graft the two estimates together,
  # preserving the middle part for grabill, and blending
  # aggressively into the young and closeout parts of Sprague
  # weights for grafting in grabill
  m                 <- nrow(pops)
  lr                <- m - 1
  fr                <- lr - 9

  # these weights do much better than linear weights.
  w10               <- exp(row(pops[1:10, , drop = FALSE])) / exp(10.1)

  # blend together young ages
  popg[1:10,]       <- w10 * popg[1:10,] + (1 - w10) * pops[1:10,]

  # blend together old ages
  popg[fr:lr,]      <- w10[10:1,] * popg[fr:lr,] + (1 - w10[10:1,]) * pops[fr:lr,]

  # ---------------------------------------------
  # now we take care of the marginal constraint problem
  # make weighting matrix
  wr                <- pops * 0 + 1
  wr[1:10,]         <- w10
  wr[fr:lr,]        <- w10[10:1,]
  wr[nrow(wr),]     <- 0

  # weighted marginal sums. The difference we need to redistribute
  redist            <- colSums(pops) - colSums(popg)

  middle.part       <- popg * wr

  # the difference to redistribute
  add.in            <- t(t(prop.table(middle.part, 2)) * redist)
  popg              <- popg + add.in
  # ---------------------------------------------
  # label dims and return
  AgeOut            <- min(Age):(min(Age) + nrow(popg) - 1)
  dim(popg)         <- NULL
  # TR: this is necessary if final age group is open,
  # but not in an age evenly divisible by 5
  names(popg)       <- min(Age):(length(popg) - 1)

  popg
}


# -------------------------------------------------------------------------- #
# Beers                                                                      #
# -------------------------------------------------------------------------- #

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
#' @description This method offers both ordinary and modified Beers splitting, with an optional \href{https://www.census.gov/data/software/dapps.html}{Demographic Analysis & Population Projection System Software} adjustment `johnson` for ages under 10.
#'
#' @inheritParams graduate
#' @param method character. Valid values are `"ord"` or `"mod"`. Default `"ord"`.
#' @param johnson  logical. Whether or not to adjust young ages according to the \href{https://www.census.gov/data/software/dapps.html}{Demographic Analysis & Population Projection System Software} method. Default `FALSE.`
#' @details Ages should refer to lower age bounds. `Value` must be labeled with ages unless `Age` is given separately. There must be at least six 5-year age groups (including the open group, 5 otherwise). If you want the `johnson` adjustment then `Value` must contain a single-year estimate of the population count in age 0. That means `Value` must come either as standard abridged or single age data.
#' 
#' `method` option `"ord"` conserves sums in 5-year age groups, whereas `"mod"` does some smoothing between 5-year age groups too, and is not constrained.
#'
#' If the highest age does not end in a 0 or 5, and `OAG == TRUE`, then the open age will be grouped down to the next highest age ending in 0 or 5. If the highest age does not end in a 0 or 5, and `OAG = FALSE`, then results extend to single ages covering the entire 5-year age group.
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
#' # graduate() deals with this automatically.
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
                           method = "ord",
                           johnson = FALSE) {

  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
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
    pop1 <- graduate_beers_johnson(Age0 = Age0,
                                   pop5 = pop5,
                                   pop1 = pop1)
  }

  # TR will this fail if Age starts with 1?
  zero              <- min(Age)
  ages              <- zero:(nrow(bm) - 1 + zero)
  names(pop1)       <- ages
  pop1
}

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
graduate_beers_johnson <- function(Age0, pop5, pop1) {

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

# -------------------------------------------------------------------------- #
# PCLM                                                                       #
# -------------------------------------------------------------------------- #

#' wrapper for \code{ungroup::pclm} method of splitting binned counts
#'
#' @description This is exactly the function \code{pclm()} from the \code{ungroup} package, except with arguments using standard \code{DemoTools} argument names.
#' @details The PCLM method can also be used to graduate rates using an offset if both numerators and denominators are available. In this case \code{Value} is the event count and \code{offset} is person years of exposure. The denominator must match the length of \code{Value} or else the length of the final single age result \code{length(min(Age):OAnew)}.  This method can be used to redistribute counts in the open age group if \code{OAnew} gives sufficient space. Likewise, it can give a rate extrapolation beyond the open age.
#' 
#' If there are 0s in `Value`, these are replaced with a small value prior to fitting. If negatives result from the pclm fit, we retry after multiplying `Value` by 10, 100, or 1000, as sometimes a temporary rescale for fitting can help performance.
#' 
#' `Age` be any age range, it does not need to start at 0.
#'
#' @inheritParams graduate
#' @param ... further arguments passed to \code{ungroup::pclm()}
#'
#' @references
#' \insertRef{pascariu2018ungroup}{DemoTools}
#' \insertRef{rizzi2015efficient}{DemoTools}
#'
#' @importFrom ungroup pclm
#' @export
#' @seealso \code{\link[ungroup]{pclm}}
#' @examples
#' a5  <- seq(0,100,by=5)
#' p5  <- pop5_mat[, 1]
#' p1  <- graduate_pclm(Value = p5, Age = a5)
#' p1s <- graduate_sprague(Value = p5, Age = a5)
#' \dontrun{
#' plot(a5, p5/5, type = "s",xlim=c(40,60),ylim=c(2000,4000))
#' lines(0:100, p1, lwd = 2, col = "red")
#' lines(0:100, p1s, lwd = 1, col = "blue",lty="8282")
#' }
#' # example of how to graduate rates by splitting deaths using population
#' # as PCLM offset
#' dth.ind <- c(49, 14, 9, 39, 60, 101, 147, 178, 177, 232)
#' pop.ind <- c(7231, 28400, 66836, 52380, 38022, 36886, 26145, 14205, 6406, 
#'              3322)
#' age <- c(0,1,5,15,25,35,45,55,65,75)
#' mx <- graduate_pclm(Value = dth.ind, 
#'               Age = age, 
#'               OAnew = 85, 
#'               offset = pop.ind)
#' \dontrun{
#'   plot(age, dth.ind / pop.ind, type = 's', log = 'y', xlim = c(0,85))
#'   lines(0:85, mx, col = "red")
#' }

graduate_pclm <- function(Value, Age, AgeInt, OAnew = max(Age), OAG = TRUE, ...) {
  
  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
  }
  
  if (OAnew > max(Age)){
    nlast    <- OAnew - max(Age) + 1
  } else {
    nlast    <- 1
  }
  a1       <- min(Age):OAnew
  DOTS     <- list(...)
  if ("offset" %in% names(DOTS)) {
    
    # offset could be one or another thing..
    lo     <- length(DOTS$offset)
    o1     <- length(a1) == lo
    o5     <- length(Value) == lo
    stopifnot(o1 | o5)
  }

  # TR 22 March 2021
  # 0s cause breakage
  # check for 0s
  ind0 <- Value == 0
  have0s <- any(ind0)
  if (have0s){
    cat("\n0s detected in Value, replacing with .01\n")
    Value[ind0] <- .01
  }
  
  A        <- pclm(x = Age, y = Value, nlast = nlast, ...)
  fac <- 1
  for (i in 1:3){
    if (any(A$fitted < 0)){
      # let's assume it's a scale issue
      fac      <- 10^i
      A        <- pclm(x = Age, y = Value * fac, nlast = nlast, ...)
    } else {
      break
    }
  }
  if (any(A$fitted < 0)){
    # TR: just let the error propagate instead of interpreting it?
    cat("\nCareful, results of PCLM produced some negatives. 
        \nWe tried rescaling inputs by as much as",fac,"\nbut alas it wasn't enough.\n")
  }
  if (fac > 1){
    cat("\nPossible small counts issue with these data and the PCLM method\nIt seems to have worked without producing negatives when fitting Value is scaled by",fac,"\nCouldn't hurt to eyeball results!\n")
  }

  B        <- A$fitted / fac
  a1.fitted <- A$bin.definition$output$breaks["left", ]
  names(B) <- a1.fitted
  # in case OAnew is lower than max(Age)
  C        <- groupOAG(Value = B, Age = a1.fitted, OAnew = OAnew)
  C
}


# -------------------------------------------------------------------------- #
# Monotonic                                                                  #
# -------------------------------------------------------------------------- #

#' Graduate age groups using a monotonic spline.
#' @description Take the cumulative sum of \code{Value} and then run a monotonic spline through it. The first differences split back single-age estimates of \code{Value}. Optionally keep the open age group untouched.
#'
#' @details The \code{"hyman"} method of \code{stats::splinefun()} is used to fit the spline because 1) it passes exactly through the points, 2) it is monotonic and therefore guarantees positive counts, and 3) it seems to be a bit less wiggly (lower average first differences of split counts). Single-age data is returned as-is. If you want to use this function as a smoother you first need to group to non-single ages. `Age` be any age range, it does not need to start at 0.
#' @inheritParams graduate
#' @return Numeric. vector of single smoothed age counts.
#' @importFrom stats splinefun
#' @references
#' \insertRef{fritsch1980monotone}{DemoTools}
#' @export
#' @examples
#' Value                <- structure(c(88623, 90842, 93439, 96325, 99281, 102051, 104351,
#'				 106555, 109170, 112188, 113582, 112614, 108904, 102622, 95867,
#'				 80874, 60196, 37523, 17927, 5642, 1110), .Names = c("0", "5",
#'				 "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60",
#'				 "65", "70", "75", "80", "85", "90", "95", "100"))
#'
#' # if the last age group is closed, then it's best to use AgeInt, otherwise,
#' # one is assumed from the age siphoned from the names attribute of Value.
#' graduate_mono(Value, OAG = FALSE)
#' # or leave open age group in tact
#' graduate_mono(Value, OAG = TRUE)
#'
#' data(pop1m_ind)
#' Value5                <- groupAges(pop1m_ind,Age=0:100,N=5) 
#' 
#' Value1 <- graduate_mono(Value = Value5, Age = names2age(Value5), OAG =TRUE)
#' 
#' \dontrun{
#'  
#'  plot(seq(0,100,5),Value5 / 5, xlab = 'Age', ylab = 'Counts', type = 's')
#'  lines(0:100,Value1)
#' }

graduate_mono   <- function(
  Value,
  Age,
  AgeInt,
  OAG = TRUE) {

  if (missing(Age) & missing(AgeInt)) {
    Age                 <- names2age(Value)
  }
  if (missing(AgeInt)) {
    # give 1 to final interval to preserve
    AgeInt              <- age2int(Age, OAG = OAG, OAvalue = 1)
  }
  if (missing(Age)) {
    Age                 <- int2age(AgeInt)
  }

  # if age is single return as-is
  if (is_single(Age)) {
    names(Value) <- Age
    return(Value)
  }

  # if last age Open, we preserve it
  if (OAG) {
    N                   <- length(Value)
    OAvalue             <- Value[N]
    Value               <- Value[-N]
    Age                 <- Age[-N]
    AgeInt              <- AgeInt[-N]
  }
  # if the final age is Open, then we should remove it and then
  # stick it back on

  AgePred               <- c(min(Age), cumsum(AgeInt) + min(Age))
  y                     <- c(0, cumsum(Value))
  AgeS                  <- min(Age):(sum(AgeInt)+ min(Age))
  # TR: changed from monoH.FC to hyman 3.3.2021
  y1                    <- splinefun(y ~ AgePred, method = "hyman")(AgeS)
  out                   <- diff(y1)
  names(out)            <- AgeS[-length(AgeS)]

  # The open age group is maintained as-is.
  if (OAG) {
    out                 <- c(out, OAvalue)
  }
  age1 <- min(Age):(min(Age) + length(out) - 1)
  names(out) <- age1
  out
}

#' blend the Sprague upper boundary age estimates into monotonic spline estimates
#'
#' @description A simple monotonic spline on the cumulative sum of population counts may return more convincing single age count estimates than the Sprague or other splitting methods. This function blends the given single age population estimates starting at \code{pivotAge}.
#'
#' @inheritParams graduate
#' @param pops optional numeric vector of single age population counts derived from \code{Value}.
#' @param pivotAge integer (default 90). Age at which to switch to spline-based estimates.
#' @param splitfun optional. The function used to create pops. Default \code{graduate_sprague}. Could also be \code{graduate_grabill}, or any other function that similarly transforms.
#' @param OAG logical (default \code{FALSE}). Would we like to re-impute the last
#' element of \code{Value} as the open age group?
#' @param ... arguments to be optionally passed to \code{splitfun()}.
#' @return numeric matrix of age by year estimates of single-age counts.
#'
#' @details The \code{pivotAge} must be at least 10 years below the maximum age detected from
#' \code{rownames(popmat)}, but not lower than 75. In the exact \code{pivotAge}, we may either take the Sprague estimates or the spline estimates, depending on which is larger, then the single-age estimates for this 5-year age group are rescaled to sum to the original total in \code{Value}. Higher ages are taken from the spline-based age splits. The spline results are derive from the \code{"hyman"} method of \code{splinefun()} on the cumulative sum of the original age grouped data. One could use this function to perform the same closeout to Grabill estimates, if these are given via the \code{pops} argument. See examples. Note that the Grabill split method mixed with this closeout will not necessarily preserve the annual totals, and this function performs to rescaling. The open age group is preserved (and must be included in \code{Value}).
#'
#' @export
#'
#' @examples
#' a5 <- as.integer(rownames(pop5_mat))
#' popvec               <- pop5_mat[,1]
#' closed.out           <- graduate_mono_closeout(Value = popvec, Age = a5, OAG = TRUE)
#' sum(closed.out) - sum(popvec)
#' graduate_mono_closeout(Value = popvec, pivotAge = 85, Age = a5, OAG = TRUE)
#' # giving a different single-age split to close out this way:
#' popg                 <- graduate_grabill(Value = popvec, Age = a5, OAG = TRUE)
#' grabill.closed.out   <- graduate_mono_closeout(Value = popvec, Age = a5, pops = popg)
#' # totals not necessarily preserved if mixed w Grabill
#' # I wouldn't recommend a rescale of the total, since the
#' # only part we mess with here is the old age section. Ergo,
#' # one may wish to instead rescale results colSums() of
#' # popg at age pivotAge and higher.
#' sum(grabill.closed.out) - sum(popvec)
#' # also works on an age-labeled vector of data

#' closed.vec           <- graduate_mono_closeout(popvec, Age = a5, OAG = TRUE)
#' # let's compare this one with sprague()
#' simple.vec           <- graduate_sprague(popvec, Age = a5, OAG = TRUE)
#' # and with a simple monotonic spline
#' mono.vec             <- graduate_mono(popvec, Age = a5, OAG = TRUE)
#' \dontrun{
#' plot(85:100,simple.vec[86:101], type = 'l',
#'  main = "In this case graduate_sprague() is the smoothest")
#' lines(85:100,closed.vec[86:101], col = "red", lwd = 2)
#' lines(85:100,mono.vec[86:101], col = "blue", lty = 2)
#' legend("topright",lty=c(1,2,2), col = c("black","red","blue"),lwd = c(1,2,1),
#' 		legend = c("graduate_sprague()","monoCloseout()", "graduate_mono()"))
#' }

graduate_mono_closeout <-
  function(Value,
           Age,
           pops,
           pivotAge = 90,
           splitfun = graduate_sprague,
           OAG = TRUE,
           ...) {
 
    names(Value)        <- Age

    if (missing(pops)) {
      pops              <- splitfun(Value, Age = Age, OAG = OAG, ...)
    }
    # get the spline population split
    AgeIn               <- Age
    # this does not smooth single ages, it only splits to single
    popmono             <- graduate_mono(Value, OAG = OAG, Age = AgeIn)

    # some age pars
    Age                 <- as.integer(names(popmono))

    # some checks on pivotAge...
    if (!(max(Age) - 10) >= pivotAge) {
      pivotAge          <- max(Age) - 10
      if (pivotAge < 75) {
        warning(
          "pivotAge wasn't in rownames(popmat), moved it to 3rd ",
          "from bottom row of popmat, but appears to be < 75 ",
          "so returning sprague() output as-is, no extra closeout performed."
        )
        return(pops)
      }
      warning("pivotAge moved to ", pivotAge, ", continued.")
    }
    # -----------------------------
    # now begin the closeout blend.
    p.i                 <- which(Age == pivotAge)
    ## substitute Sprague interpolation if > pchip for better blending of the two series
    pop.c               <- popmono[p.i:(p.i + 4)]
    ind                 <- pops[p.i] > pop.c[1]
    if (ind){
      pop.c[1]       <- pops[p.i]
    }


    ## adjust back on initial pop 5x5 for age 90-94
    ## proportional distribution
    pop.c[is.na(pop.c)] <- 0
    prop                <- rescale_vector(pop.c, scale = 1)
    pivot5              <- Value[which(AgeIn == pivotAge)]
    pop.c               <- prop * pivot5
    ## append the remaining of the age groups (except last open age)
    ## 95-99 onward
    m                   <- length(pops)
    pop.c               <- c(pop.c, popmono[(p.i + 5):m])
    ## append Sprague interpolation before age 90
    pop.c               <- c(pops[1:(p.i - 1)], pop.c)

    ## deal with negative values if applicable (but in principle should not be happening)
    pop.c[pop.c < 0]    <- 0

    # label and return
    #dimnames(pop.c)    <- list(Age1, colnames(popmat))
    names(pop.c)     <- Age

    pop.c
  }

#' Graduate grouped data
#'
#' @description A wrapper function for several graduation methods, primarily for count data (\code{"sprague"}, \code{"beers(ord)"}, \code{"beers(mod)"}, \code{"grabill"}, \code{"mono"} (Monotonic spline), \code{"uniform"}, \code{"pclm"}), but also with one (\code{"pclm"}) with an option for graduating rates if both event counts and population at risk are available.
#'
#' @details \code{"sprague"}, \code{"beers(ord)"}, \code{"beers(mod)"} methods require original data to be in uniform five-year age groups. If they are not (for example, the infant group is separate) then they are grouped to uniform width prior to splitting. If you want to keep the original infant count in output, then specify \code{keep0 = TRUE}. In this case, it is imputed, and ages 1-4 are rescaled, which may introduce a discontinuity in results from age 4 to 5. \code{keep0 = TRUE} may also be desired along with \code{method = "pclm"}.
#'
#' Some methods are constrained, others not, and others are optionally constrained. If this is required, then this function can be followed up with \code{rescaleAgeGroups()}, which may have the effect of breaking continuity in smooth results. This is inconsequential for downstream demography, but if this aesthetic side effect is undesired, then try one of the constrained methods: \code{"sprague"}, \code{"mono"}, \code{"pclm"} (with \code{control = list(lambda = 1/1e7) specified} or similar).
#'
#' Beers may either be ordinary \code{"beers(ord)"} or modified \code{"beers(mod)"}, and either can pass on the optional argument \code{johnson = TRUE} if desired (this has a different distribution pattern for young ages, \code{FALSE} by default). If \code{method = "beers"} is given, then \code{"beers(ord)"} is used.
#'
#' This wrapper standardizes some inconsistencies in how open ages are dealt with. For example, with the \code{"pclm"} method, the last age group can be redistributed over a specified interval implied by increase \code{OAnew} beyond the range of \code{Age}. To get this same behavior from \code{"mono"}, or \code{"uniform"} specify \code{OAG = FALSE} along with an appropriately high \code{OAnew} (or integer final value of \code{AgeInt}.
#'
#' \code{OAnew} cannot be higher than \code{max(Age)+4} for \code{"sprague"} or \code{"beers"} methods. For \code{"uniform","mono","pclm"} it can be higher than this, and in each case the open age group is completely redistributed within this range, meaning it's not really open anymore.
#'
#' For all methods, negative values are detected in output. If present, we deal with these in the following way: we take the geometric mean between the given output (with negative imputed with 0s) and the output of \code{graduate_mono()}, which is guaranteed non-negative. This only affects age groups where negatives were produced in the first pass. In our experience this only arises when using Sprague, Beers, or Grabill methods, whereas all others are guaranteed non-negative.
#'
#' For any case where input data are in single ages, constraining results to sum to values in the original age groups will simply return the original input data, which is clearly not your intent. This might arise when using graduation as an implicit two-step smoother (group + graduate). In this case, separate the steps, first group using \code{groupAges()} then use \code{graduate(..., constrain = TRUE)}.
#'
#' @param Value numeric vector, presumably counts in grouped ages
#' @param Age integer vector, lower bounds of age groups
#' @param AgeInt integer vector, age interval widths
#' @param OAG logical, default = \code{TRUE} is the final age group open?
#' @param OAnew integer, optional new open age, higher than \code{max(Age)}. See details.
#' @param method character, either \code{"sprague"}, \code{"beers(ord)")}, \code{"beers(mod)")}, \code{"mono")}, \code{"uniform")}, or \code{"pclm"}
#' @param keep0 logical. Default \code{FALSE}. If available, should the value in the infant age group be maintained, and ages 1-4 constrained?
#' @param constrain logical. Default \code{FALSE}. Should output be constrained to sum within the input age groups?
#' @param ... extra arguments passed to \code{graduate_beers()} or \code{graduate_pclm()}
#' @seealso \code{\link{graduate_sprague}}, \code{\link{graduate_beers}}, \code{\link{graduate_uniform}}, \code{\link{graduate_mono}}, \code{\link{graduate_pclm}}, \code{\link{graduate_grabill}}
#' @export
#' @references
#' \insertRef{pascariu2018ungroup}{DemoTools}
#' \insertRef{rizzi2015efficient}{DemoTools}
#' \insertRef{sprague1880explanation}{DemoTools}
#' \insertRef{shryock1973methods}{DemoTools}
#' \insertRef{siegel2004methods}{DemoTools}
#' \insertRef{beers1945modified}{DemoTools}
#'
#' @examples
#' Value <- pop5_mat[, 1]
#' Value <- c(10000,44170,Value[-1])
#' Age   <- sort(c(1,seq(0,100,by=5)))
#'
#' graduate(Value, Age, method = "sprague")
#' graduate(Value, Age, method = "sprague", keep0=FALSE)
#'
#' graduate(Value, Age, method = "beers(ord)")
#' graduate(Value, Age, method = "beers(ord)", keep0=TRUE)
#' graduate(Value, Age, method = "beers(ord)", keep0=TRUE, johnson = TRUE)
#'
#' graduate(Value, Age, method = "beers(mod)")
#' graduate(Value, Age, method = "beers(mod)", keep0=TRUE)
#' graduate(Value, Age, method = "beers(mod)", keep0=TRUE, johnson = TRUE)
#'
#' graduate(Value, Age, method = "mono")
#' graduate(Value, Age, method = "mono", keep0=TRUE)
#'
#' graduate(Value, Age, method = "uniform")
#'
#' graduate(Value, Age, method = "pclm")
#' graduate(Value, Age, method = "pclm", keep0=TRUE)
#' # pclm can also graduate rates if both
#' # numerators and denominators are on hand:
#' Exposures <- c(100958,466275,624134,559559,446736,370653,301862,249409,
#'                247473,223014,172260,149338,127242,105715,79614,53660,
#'                31021,16805,8000,4000,2000,1000)
#'
#' Deaths <- c(8674,1592,618,411,755,1098,1100,1357,
#'             1335,3257,2200,4023,2167,4578,2956,4212,
#'             2887,2351,1500,900,500,300)
#' Age    <- c(0, 1, seq(5, 100, by = 5))
#' AgeInt <- c(diff(Age), NA)
#'
#' # exclude infants for better fit.
#' mx    <- graduate(
#'            Value = Deaths[-1], Age = Age[-1],
#'            AgeInt = AgeInt[-1], OAG = TRUE,
#'            OAnew = 110, offset = Exposures[-1],
#'            method = "pclm")
#' mx_sm <- graduate(
#'            Value = Deaths[-1], Age = Age[-1],
#'            AgeInt = AgeInt[-1], OAG = TRUE,
#'            OAnew = 110, offset = Exposures[-1],
#'            method = "pclm", control = list(lambda = 1e7))
#'
#' \dontrun{
#' plot(Age,
#'      Deaths / Exposures,
#'      type = 's', log = 'y',
#'      main = "Underlying data have differential heaping on 0s and 5s")
#' lines(1:110, mx)
#' lines(1:110, mx_sm, col = "blue")
#' legend("bottomright",
#'        col = c("black","blue"),
#'        lty = c(1, 1),
#'        legend = c("lambda optimized (almost constrained)",
#'                   "higher lambda = smoother")
#'        )
#'   }

graduate <- function(Value,
                     Age,
                     AgeInt = age2int(Age),
                     OAG = TRUE,
                     OAnew = max(Age),
                     method = c("sprague",
                                "beers(ord)",
                                "beers(mod)",
                                "grabill",
                                "pclm",
                                "mono",
                                "uniform"),
                     keep0 = FALSE,
                     constrain = FALSE,
                     ...) {
  method <- tolower(method)
  if (method == "beers") {
    method == "beers(ord)"
  }
  # validate method choice
  method <- match.arg(method)

  # TR: those methods that require AgeInt may depend on final value not being NA
  # even if it's an open age, in essence, keep all value inside this same "single"
  # age. NA coding isn't the best choice here, but we anticipate this and assign
  # 1 to AgeInt[length(AgeInt)] IFF OAG & max(Age) == OAnew
  # This is inconsequential for those methods that don't use AgeInt
  N <- length(AgeInt)
  if (OAG & is.na(AgeInt[N])) {
    nlast     <- OAnew - max(Age) + 1
    AgeInt[N] <- nlast
  }

  # Sprague in strict 5-year age groups
  if (method == "sprague") {
    out <- graduate_sprague(Value, Age = Age, OAG = OAG)
  }

  # Grabill in strict 5-year age groups
  if (method == "grabill") {
    out <- graduate_grabill(Value, Age = Age, OAG = OAG)
  }
  # Beers in strict 5-year age groups
  if (grepl("beers", method)) {
    if (grepl("ord", method)) {
      out <- graduate_beers(Value,
                   Age = Age,
                   OAG = OAG,
                   method = "ord",
                   ...)
    }
    if (grepl("mod", method)) {
      out <- graduate_beers(Value,
                   Age = Age,
                   OAG = OAG,
                   method = "mod",
                   ...)
    }
  }

  # inconsistent ways of dealing with top age group.
  # if !OAG, then there is potential for inconsistency
  # between (min(Age) + sum(AgeInt)) and OAnew,
  # namely, OAnew can increase or decrease this age IFF
  # OAG, but if there is an inconsistency then we need
  # to either throw an error or declare a preference.
  # TR: preference should go to AgeInt, because it refers
  # to data observations, whereas OAnew refers to a desired
  # change. So actually no change is needed?


  # Uniform respects irregular intervals
  if (method == "uniform") {
    OAvalue <- OAnew - max(Age) + 1
    out <- graduate_uniform(
      Value = Value,
      Age = Age,
      AgeInt = AgeInt,
      OAG = OAG,
      OAvalue = OAvalue
    ) # OAvalue only if OAG and
    # extrapolation desired?
  }

  # Mono respects irregular intervals
  if (method == "mono") {
    out <- graduate_mono(
      Value = Value,
      Age = Age,
      AgeInt = AgeInt,
      OAG = OAG
    ) # doesn't allow for extension?
    # actually it does IFF !OAG & !is.na(AgeInt[length(AgeInt)])
  }

  # PCLM respects irregular intervals
  if (method == "pclm") {
    out <- graduate_pclm(Value = Value,
                         Age = Age,
                         OAnew = OAnew,
                         ...)
  }

  # handle infant age for sprague and beers, if required.
  if (keep0) {
    if (Age[1] == 0 & AgeInt[1] == 1) {
      V0       <- Value[1]
      V5       <- sum(Value[Age < 5])
      V4       <- V5 - V0
      a1       <- names2age(out)
      ind      <- a1 < 5 & a1 > 0
      out[ind] <- rescale_vector(out[ind], scale = V4)
      out[1]   <- V0
    }
  }

  n  <- length(out)
  a1 <- min(Age):(min(Age) + n - 1)

  # detect negatives. Have default option to replace.
  # Favor quick over perfect, since this only can arise
  # in Sprague, Beers, or Grabill, which are quick. In this
  # case the age group with negatives will no longer be constrained
  # to sum to the same original value.
  ind0 <- out < 0
  if (any(ind0)){
    # which
  
    agen         <- rep(Age, times = AgeInt)
    problem.ages <- agen[ind0]
    out[ind0]    <- 0
    # well, this is maybe to complicated, can swap w
    # different trick
    outm <- graduate_mono(
              Value = Value,
              Age = Age,
              AgeInt = AgeInt,
              OAG = OAG)

    swap <- agen %in% problem.ages

    out.swap <- interp(cbind(out[swap],outm[swap]),
                       datesIn = c(10,20),
                       datesOut = 15,
                       method = "power",
                       power = 1/2)

    out[swap] <- out.swap
    dim(out)  <- NULL
  }

  # option to contrain to sum to original age groups
  if (constrain){
    out <- rescaleAgeGroups(Value1 = out,
                            AgeInt1 = rep(1, n),
                            Value2 = Value,
                            AgeInt2 = AgeInt,
                            splitfun = graduate_uniform,
                            recursive = FALSE)
    out[is.nan(out)] <- 0
  }

  # last min names assure
  names(out) <- a1

  out
}
