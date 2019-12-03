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
#' @details  Dimension labelling is necessary. There must be at least six age groups (including the open group). One year of data will work as well, as long as it's given as a single-column matrix. Data may be given in either single or grouped ages. If the highest age does not end in a 0 or 5, and \code{OAG == TRUE}, then the open age will be grouped down to the next highest age ending in 0 or 5. If the highest age does not end in a 0 or 5, and \code{OAG == FALSE}, then results extend to single ages covering the entire 5-year age group.
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
  OAG = TRUE) {
  
  if (as.character(match.call()[[1]]) == "grabill") {
    warning("please use graduate_grabill() instead of grabill().", call. = FALSE)
  }
  
  punif1       <- graduate_uniform(
                    Value = Value, 
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
  names(popg)       <- a1
  
  popg
}

#' @export
#' @rdname graduate_grabill
grabill <- graduate_grabill