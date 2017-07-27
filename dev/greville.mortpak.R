lt <- readRDS("C:/Users/pgerl/Dropbox/bayesPop/bayesPop2017/LTWorldMale.rds")

greville.mortpak <- function(QxMx, inputQxMx, Sex) {
  # QxMx is either Qx or Mx vector of values  
  # with inputQxMx = 1 for Qx as input for QxMx, and 2 for Mx otherwise like in Mortpak LIFETB
  # with Sex =1 for Male and 2 for Female
 if (Sex==1){ # Male under age 5 based on CD-West
     if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
       if (QxMx[1] < 0.1) {
         ax0 <- 0.0425 + 2.875 * QxMx[1]
         ax1 <- 1.653 - 3.013 * QxMx[1]
       } else 
       {
         ax0 <- 0.330
         ax1 <- 1.352
       }
     }
     if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
       if (QxMx[1] < 0.107) {
           ax0 <- 0.045 + 2.684 * QxMx[1]
           ax1 <- 1.651 - 2.816 * QxMx[1]
       } else 
         {
           ax0 <- 0.330
           ax1 <- 1.352
         }
     }
  }
 if (Sex==2){ # Female under age 5 based on CD-West
    if (inputQxMx == 1) {  #if input is Qx (Table on p. 20 in Coale and Demeny 1983 for CD-West)
      if (QxMx[1] < 0.1) {
        ax0 <- 0.050 + 3.000 * QxMx[1]
        ax1 <- 1.524 - 1.625 * QxMx[1]
      } else 
      {
        ax0 <- 0.350
        ax1 <- 1.361
      }
    }
    if (inputQxMx == 2) {  #if input is Mx (Table 3.3, p. 48 in Preston et al 2001)
      if (QxMx[1] < 0.107) {
        ax0 <- 0.053 + 2.800 * QxMx[1]
        ax1 <- 1.522 - 1.518 * QxMx[1]
      } else 
      {
        ax0 <- 0.350
        ax1 <- 1.361
      }
    }
  }
  
  ## Greville (based on Mortpak LIFETB) for other ages
  ax <- c(NA, 2.5 - (25 / 12) * (QxMx[-c(1, length(QxMx))] - 0.1 * log(QxMx[3:length(QxMx)] / QxMx[1:(length(QxMx)-2)])),
            1/QxMx[length(QxMx)])
  ax[1:4] <- c(ax0, ax1, 2.5, 2.5)
  ax
}

QxMx <- qx

n <- c(1, 4, rep(5, length(QxMx)-2))

mx.from.qx <- function(qx, ax, n) {
    qx / (n - (n - ax) * qx)      }

qx.from.mx <- function(mx, ax, n)    {
    (n * mx) / (1.0 + ((n - ax) * mx)) }

qx <- lt$qx

## initialize ax using Qx
ax1 <- greville.mortpak(qx, 1, 1)
axi <- ax1

## iterate until output qx match input qx using mxi with axi and qx as input
smsq <- 99999
epsilon <- 1.00E-16
while(smsq > epsilon) {
  mxi   <- mx.from.qx(qx, axi, n)
  axi   <- greville.mortpak(mxi, 2)
  qxnew <- qx.from.mx(mxi, axi, n)
  smsq  <- sum((qxnew - qx)^2)
  print(smsq)
}

mx <- mxi
ax <- axi
## note the last mx for the open age group needs to be computed differently, or rather mx only up to age 95-99 used to fit Kannisto, and the open age group discarded.
mx[length(mx)] <- NA
ax[length(mx)] <- NA

check <- data.frame(cbind(mx, qxnew, qx, ax))

#install.packages("bayesPop")




