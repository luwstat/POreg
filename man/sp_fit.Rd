\name{sp_fit}
\alias{sp_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression analysis of arbitrarily censored data under the proportional odds model
}
\description{
Fits the semi-parametric maximum likelihood estimator for univariate arbitrarily censored data using EM algorithm.
}
\usage{
sp_fit <- function(L,R,x,order = 3,equal_space = T, nknot,myknots, conv_cri = 1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{Lower bound of the observed intervals, a length n vector.}
  \item{R}{Upper bound of the observed intervals, a length n vector.}
  \item{x}{Covariate matrix, an n x p matrix.}
  \item{order}{The order of I-spline basis function.}
  \item{equal_space}{Logical, if TRUE knots are placed evenly across the observed intervals based on the input data set; if FALSE knots should be specified by the user. see myknots. The default is TRUE.}
  \item{nknot}{The total number of knots.}
  \item{myknots}{Knots specified by the user.}
  \item{conv_cri}{The convergence criterion of the EM algorithm.}
}
\details{
L and R are left and right ends of the censoring interval. This allows for left censored (L == 0), right censored (R == inf), uncensored (L == R) along with interval censored observations (0 < L < R < inf).

The EM algorithm converges when the sum of square of the difference in the parameter estimates is less than conv_cri.

}
\value{
\item{beta}{Estimates of regression coefficients.}
\item{spline_coef}{Estimates of spline coefficients.}
\item{beta_se}{The standard errors of regression coefficients estimate.}
\item{CI}{Confidence intervals for the regression coefficients.}
}
\references{
Lu Wang & Lianming Wang. Regression analysis of arbitrarily censored data under the proportional odds models
}
\author{
Wang,L. and Wang,L.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(icenReg)
data("IR_diabetes")
gender <- as.matrix(as.numeric(IR_diabetes[,'gender'] == 'male'))
sp_fit(L = IR_diabetes[,'left'],R = IR_diabetes[,'right'],x = gender,order = 3,nknot = 10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
