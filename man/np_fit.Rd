\name{np_fit}
\alias{np_fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-parametric estimator for arbitrarily censored data under the proportional odds model
}
\description{
Fits the non-parametric maximum likelihood estimator (NPMLE) for univariate arbitrarily censored data using I-spline.
}
\usage{
np_fit(L,R,order = 3,equal_space = T,nknot, myknots,conv_cri = 1e-9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{Lower bound of the observed intervals, a length n vector.}
  \item{R}{Upper bound of the observed intervals, a length n vector.}
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
\item{spline_coef}{Estimates of spline coefficients.}
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
data("marijuana")

np_fit(marijuana[,1],marijuana[,2],nknot = 10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }% use one of  RShowDoc("KEYWORDS")
\keyword{ }% __ONLY ONE__ keyword per line
