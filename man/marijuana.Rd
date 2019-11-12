\name{marijuana}
\alias{marijuana}
\docType{data}
\title{
First Use of Marijuana
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("marijuana")}
\format{
  A data set with 191 observations of the age when the teenager first used marijuana. It contains the following 2 variables.
  \describe{
    \item{\code{L}}{a numeric vector of the lower bound of the observed interval}
    \item{\code{R}}{a numeric vector of the upper bound of the observed interval}
  }
}
\details{
In 1975, on the Stanford-Palo Alto Peer Counseling Program, Hamburg et al. studied drug use in a representative sample of suburban junior and senior high school students. In 1987, Turnbull & Weiss extracted a set of failure time data of marijuana use from Hamburg et al.'s study which is the data set here.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Turnbull, B. W. & Weiss, L. (1978), A likelihood ratio statistic for testing goodness of fit with randomly censored data, Biometrics 34, 367-375.

Hamburg, A. B., Kraemer, C. & Jahnke, W. (1975), A hierarchy of drug use in adolescence: Behavioral and attitudinal correlates of substantial drug use, Am J Psychiatry 132, 1155-1163.

Sun, J. (2006), The statistical analysis of interval-censored data, Springer.
}
\examples{
data(marijuana)

fitted <- np_fit(marijuana[,1],marijuana[,2],nknot = 10,order = 3)

#time points
t <- seq(0,max(setdiff(c(marijuana[,1],marijuana[,2]),Inf)) + 0.01,0.05)
#calculate baseline odds
g <- fitted$spline_coef
knots <- fitted$knots
b <- Ispline(t,3,knots)
yy <- t(as.matrix(g)) %*% b
plot(t,1/(1+yy),type = "l", lty = 3,col = 'blue', ylab = "S(t)")
}
\keyword{datasets}
