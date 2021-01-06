\name{utility}
\docType{data}
\alias{utility}

\title{Fossil Fuel Fired Steam Electric Power Plants}

\description{
   The data frame \code{utility} is a panel data set 
   of 72 investor-owned fossil fuel fired steam electric power plants
   in the United States over eleven years (1986-1996).
   This data set has 791 observations and, thus, is almost balanced
   (only one firm-year observation is missing).
   This data set is used as an example in 
   Kumbhakar, Wang, and Horncastle (2015).
}

\usage{data(utility)}

\format{
   This data frame contains the following columns/variables:
   \describe{
      \item{firm}{firm number (1-72).}
      \item{year}{last two digits of the year (86-96).}
      \item{y}{net steam electric power generation in MWh.}
      \item{regu}{a dummy variable indicating whether the firm is regulated.}
      \item{k}{capital input quantity (estimate of capital cost).}
      \item{labor}{quantity index of labor and maintenance input
         (costs of labor and maintenance divided by \code{wl}).}
      \item{fuel}{quantity index of fuel input.
         (costs of fuel divided by \code{wf}).}
      \item{wl}{cost-share weighted price of labor and maintenance.}
      \item{wf}{average price of fuel (coal, oil and gas) in USD per BTU.}
      \item{wk}{price of the capital input.}
   }
}

\details{
   This data set is a revised version of the data set
   that is used in Kumbhakar and Wang (2006).
   A detailed description of the data set and the variables
   can be found in Section 4 of Kumbhakar and Wang (2006).
}

\source{
   Companion website to Kumbhakar, Wang, and Horncastle (2015):
   \url{https://sites.google.com/site/sfbook2014/}.
}

\references{
   Kumbhakar, Subal C. and Hung-Jen Wang (2006):
   Estimation of Technical and Allocative Inefficiency: A Primal System Approach,
   Journal of Econometrics 134(2), p. 419-440.
   (\doi{10.1016/j.jeconom.2005.07.001})

   Kumbhakar, Subal C., Hung-Jen Wang, and Alan P. Horncastle (2015): 
   A Practitioner's Guide to Stochastic Frontier Analysis Using Stata. 
   Cambridge University Press. 
   (\doi{10.1017/CBO9781139342070})
}

\keyword{datasets}
