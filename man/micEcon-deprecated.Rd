\name{micEcon-deprecated}
\alias{writeFront41in}
\alias{readFront41out}
\title{Deprecated Functions in the micEcon package}

\description{
   These functions are provided for compatibility with older versions
   of the micEcon package only,
   and may be defunct as soon as the next release.
}

\usage{
   writeFront41in( data, crossSectionName, timePeriodName,
   yName, xNames = NULL, zNames = NULL,
   translog = FALSE, quadHalf = TRUE,
   functionType = 1, modelType = 1, logDepVar = TRUE, mu = FALSE, eta = FALSE,
   insFile = "front41.ins", dtaFile = sub( "\\.ins$", ".dta", insFile ),
   outFile = sub( "\\.ins$", ".out", insFile ) )

   readFront41out( file = "front41.out", translog = FALSE )
}

\arguments{
   \item{data}{data frame that contains the data.}
   \item{crossSectionName}{string: name of the cross section identifier.}
   \item{timePeriodName}{string: name of the time period identifier.}
   \item{yName}{string: name of the endogenous variable.}
   \item{xNames}{a vector of strings containing the names of the X variables
      (exogenous variables of the production or cost function).}
   \item{zNames}{a vector of strings containing the names of the Z variables
      (variables explaining the efficiency level).}
   \item{translog}{logical. Should a translog function be estimated?
      / Is the estimated function a translog function?
      (not implemented in \code{readFront41out} yet).}
   \item{quadHalf}{logical. Multiply the quadratic terms of the translog
      function by one half?}
   \item{functionType}{function type: either 1 for 'production function' or
      2 for 'cost function'.}
   \item{modelType}{model type: either 1 for 'Error Components Frontier' or
      2 for 'Tech. Eff. Effects Frontier'.}
   \item{logDepVar}{logical. Is the dependent variable logged.}
   \item{mu}{logical. Should a 'mu' (if modelType = 1) or a delta0
      (if modelType = 2) be included in the estimation.}
   \item{eta}{logical. Should an 'eta' be included in the estimation
      (only if modelType = 1).}
   \item{insFile}{name of the instruction file.}
   \item{dtaFile}{name of the data file.}
   \item{outFile}{name of the output file.}
   \item{file}{character variable with the name of the file to read.}
}

\details{
   \code{writeFront41in} is superseded by \code{front41WriteInput} and
   \code{readFront41Out} is superseded by \code{front41ReadOutput}.
}

\value{
   \code{writeFront41in} writes an instruction file and an data file to disk.

   \code{readFront41out} returns a list of class \code{front41Output}
   containing following objects:
   \item{version}{the version of Frontier 4.1 that produced the output.}
   \item{insFile}{name of the instruction file used by Frontier 4.1.}
   \item{dtaFile}{name of the data file used by Frontier 4.1.}
   \item{modelType}{model type: either 1 for 'Error Components Frontier' or
      2 for 'Tech. Eff. Effects Frontier'.}
   \item{modelTypeName}{model type: 'Error Components Frontier' or 'Tech. Eff.
      Effects Frontier'.}
   \item{functionType}{function type: either 1 for 'production function' or
      2 for 'cost function'.}
   \item{functionTypeName}{function type: 'production function' or
      'cost function'.}
   \item{logDepVar}{logical. Is the dependent variable logged.}
   \item{olsResults}{results of the OLS estimation.}
   \item{nXvars}{number X variables (exogenous variables of the production
      or cost function.}
   \item{olsLogl}{log likelihood value of the OLS estimation.}
   \item{gridResults}{results of the grid search.}
   \item{mleResults}{results of the maximum likelihood estimation.}
   \item{mleLogl}{log likelihood value of the maximum likelihood estimation.}
   \item{mleCov}{coefficient covariance matrix of the maximum likelihood
      estimation.}
   \item{efficiency}{technical efficiency estimates.}
}

\references{
   Coelli, T. (1996) A Guide to FRONTIER Version 4.1: A Computer
      Program for Stochastic Frontier Production and Cost Function
      Estimation, CEPA Working Paper 96/08,
      \url{http://www.uq.edu.au/economics/cepa/frontier.htm},
      University of New England.
}

\seealso{\code{\link{front41WriteInput}}, \code{\link{front41ReadOutput}}}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\keyword{models}

