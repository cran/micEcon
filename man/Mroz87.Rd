\name{Mroz87}
\docType{data}
\alias{Mroz87}

\title{U.S. Women's Labor Force Participation}

\description{
   The \code{Mroz87} data frame contains data about 753 married women.
   These data are collected within the "Panel Study of Income Dynamics" (PSID).
   Of the 753 observations, the first 428 are for women with positive hours
   worked in 1975, while the remaining 325 observations are for women who
   did not work for pay in 1975.  A more complete discussion of the data is
   found in Mroz (1987), Appendix 1.
}

\usage{data(Mroz87)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{lfp}{Dummy variable for labor-force participation.}
      \item{hours}{Wife's hours of work in 1975.}
      \item{kids5}{Number of children 5 years old or younger.}
      \item{kids618}{Number of children 6 to 18 years old.}
      \item{age}{Wife's age.}
      \item{educ}{Wife's educational attainment, in years.}
      \item{wage}{Wife's average hourly earnings, in 1975 dollars.}
      \item{repwage}{Wife's wage reported at the time of the 1976 interview.}
      \item{hushrs}{Husband's hours worked in 1975.}
      \item{husage}{Husband's age.}
      \item{huseduc}{Husband's educational attainment, in years.}
      \item{huswage}{Husband's wage, in 1975 dollars.}
      \item{faminc}{Family income, in 1975 dollars.}
      \item{mtr}{Marginal tax rate facing the wife.}
      \item{motheduc}{Wife's mother's educational attainment, in years.}
      \item{fatheduc}{Wife's father's educational attainment, in years.}
      \item{unem}{Unemployment rate in county of residence, in percentage points.}
      \item{city}{Dummy variable = 1 if live in large city, else 0.}
      \item{exper}{Actual years of wife's previous labor market experience.}
      \item{nwifeinc}{Non-wife income.}
      \item{wifecoll}{Dummy variable for wife's college attendance.}
      \item{huscoll}{Dummy variable for husband's college attendance.}
   }
}

\source{
   Mroz, T. A. (1987)
   The sensitivity of an empirical model of married women's hours of work to
   economic and statistical assumptions.
   \emph{Econometrica} \bold{55}, 765--799.

   PSID Staff, The Panel Study of Income Dynamics,
   Institute for Social ResearchPanel Study of Income Dynamics,
   University of Michigan, \url{http://psidonline.isr.umich.edu}.
}

\references{
   Fox, J. (2004)
   \code{car}: The Companion to Applied Regression Libary -
   An R and S-PLUS Companion to Applied Regression,
   \url{http://socserv.mcmaster.ca/jfox/Books/Companion/car.html}.
}

\keyword{datasets}
