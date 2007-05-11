\name{readFront41out}
\alias{readFront41out}
\title{Read output of Frontier 4.1}

\description{
   Read the output file of Coelli's program Frontier 4.1 that
   performs stochastic frontier analysis.
}

\usage{
   readFront41out( file = "front41.out", translog = FALSE )
}

\arguments{
   \item{file}{character variable with the name of the file to read.}
   \item{translog}{logical. Is the estimated function a translog function?
      (not implemented yet).}
}

\value{
   a list of class \code{front41out} containing following objects:
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
   Battese, G.E. and T. Coelli (1992), Frontier production functions,
      technical efficiency and panel data: with application to paddy
      farmers in India. \emph{Journal of Productivity Analysis}, 3, 153-169.

   Battese, G.E. and T. Coelli (1995), A model for technical inefficiency effects
      in a stochastic frontier production function for panel data.
      \emph{Empirical Economics}, 20, 325-332.

   Coelli, T. (1996) A Guide to FRONTIER Version 4.1: A Computer
      Program for Stochastic Frontier Production and Cost Function
      Estimation, CEPA Working Paper 96/08,
      \url{http://www.uq.edu.au/economics/cepa/frontier.htm},
      University of New England.
}

\seealso{\code{\link{writeFront41in}}}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # read the output file that is provided with Frontier 4.1
   outFile <- system.file( "front41/EG1.OUT", package = "micEcon" )
   sfa <- readFront41out( outFile )
   sfa$mleResults
   sfa$efficiency

   # perform an SFA and read the output
   data( Coelli )
   Coelli$logOutput  <- log( Coelli$output )
   Coelli$logCapital <- log( Coelli$capital )
   Coelli$logLabour  <- log( Coelli$labour )

   writeFront41in( Coelli, "firm", "time", "logOutput",
      c( "logCapital", "logLabour" ) )

   \dontrun{
   system( "wine front41.exe" )
   sfa <- readFront41out()
   sfa$mleResults
   sfa$efficiency
   }
}

\keyword{models}

