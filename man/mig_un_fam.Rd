% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mig_un_fam.R
\name{mig_un_fam}
\alias{mig_un_fam}
\title{Net migration by age for an UN family}
\usage{
mig_un_fam(NM, family, Single = TRUE, OAnew = 100)
}
\arguments{
\item{NM}{numeric. Total net migration to distribute between ages and sex.}

\item{family}{character. Could be "Family", "Female Labor", "Male Labor".}

\item{Single}{logical. Results by simple age. Default \code{TRUE}.
Typically from pre-working age and working age parts of in Roger-Castro formula.}

\item{OAnew}{The age from which to group all ages into an open ended age group.
By default it is set to 100, so it groups all ages up to 120, which is the
maximum age.}
}
\value{
List with
\itemize{
\item \code{params_RC} \code{data.frame}. Roger-Castro parameters in a \code{data.frame}. Same as \code{mig_un_params} data.
\item \code{net_migr} \code{data.frame}. Net migrants by age and sex for the chosen family.
}
}
\description{
Given a total net migration,
calculate the net migration age schedule based on the Rogers and Castro formula for UN families.
}
\examples{
# 10000 net migrants, comparing two possible families
nm1 <- mig_un_fam(NM = 10000, family = "Male Labor", OAnew = 100)
nm2 <- mig_un_fam(NM = 10000, family = "Family", OAnew = 100)
# See the female profile in for these models:
\dontrun{
plot(nm1$net_migr$age[nm1$net_migr$sex=="Female"],
     nm1$net_migr$nm[nm1$net_migr$sex=="Female"],
     xlab="Age",ylab="nm",ylim=c(0,300))
points(nm2$net_migr$age[nm2$net_migr$sex=="Female"],
       nm2$net_migr$nm[nm2$net_migr$sex=="Female"], col=2)
}
}
