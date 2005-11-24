\name{probit}
\alias{probit}
\title{Probit binary choice model}
\description{
  The functions calculates probit binary choice model
}
\usage{
probit(formula, b0 = NULL, data = sys.frame(sys.parent()),
       x=FALSE, y = FALSE, model = FALSE, method="ML", ... )
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit in the
    form \code{response ~ explanatory variables}.  The
    response variable should be either logical or numeric with values
    only 0-1.  See
    e.g. \code{\link{lm}} for details.}
  \item{b0}{inital value of parameters.}
  \item{data}{an optional data frame containing the variables in the
    model. If not found in data, the
    variables are taken from environment(formula), typically the
    environment from which \code{probit} is called.}
  \item{x, y, model}{logicals. If TRUE the corresponding components of
    the fit (the
    model matrix, the response, the model frame) are returned.}
  \item{method}{the method to be used; for fitting,
    currently only method = "ML" (Maximum Likelihood) is supported;
    method = "model.frame" returns the model frame (the same as
    with model = TRUE, see below).}
  \item{\dots}{further arguments for \code{maxLik}.}
}
\details{
  The model is estimated using Maximum Likelihood method.
}
\value{
  An object of class "probit".  It is a list with following components:
  \item{LRT}{Likelihood ration test.  The full model is tested against
    H0: the parameters have no effect on the result.  This is a list
    with components
    \item{LRT}{The LRT value}
    \item{df}{Degrees of freedom for LRT (= df of the model - 1)}
    LRT is distributed by chi2(df).}
  \item{NParam}{Number of parameters of the model including constant}
  \item{NObs}{Number of the observations}
  \item{N1}{Number of observations with non-zero (true) response}
  \item{N0}{Number of observations with zero (false) response}
  \item{df}{Number of free parameters}
  \item{x}{if requested, the model matrix used.}
  \item{y}{if requested, the response used.}
  \item{model}{if requested, the model frame used.}
  Other components are inherited from \code{\link{maxLik}}.
}
\references{W. Greene: Advanced Econometrics, 2000, p ....}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{maxLik}} for ready-packaged likelihood maximisation
  routines and methods, \code{\link{glm}} for generalised linear models,
  including probit, \code{\link{binomial}}.
}
\examples{
## A simple MC trial: note probit assumes normal errors
x <- runif(100)
e <- 0.5*rnorm(100)
y <- x + e
summary(probit((y > 0) ~ x))
## female labour force participation probability
data(Mroz87)
Mroz87$kids <- Mroz87$kids5 > 0 | Mroz87$kids618 > 0
Mroz87$age30.39 <- Mroz87$age < 40
Mroz87$age50.60 <- Mroz87$age >= 50
summary(probit(lfp ~ kids + age30.39 + age50.60 + educ + hushrs +
               huseduc + huswage + mtr + motheduc, data=Mroz87))
}
\keyword{models}
\keyword{regression}
\keyword{nonlinear}
