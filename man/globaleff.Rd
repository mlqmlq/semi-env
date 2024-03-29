% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Globaleff.R
\name{globaleff}
\alias{globaleff}
\title{Estimate enveloped central space via global efficient estimator}
\usage{
globaleff(X, Y, u, gamma_init = "notgiven")
}
\arguments{
\item{X}{predictor matrix}

\item{Y}{Response matrix}

\item{u}{Dimension of Envelope of selected}

\item{gamma_init}{Initial value of gamma, defaultly obtaining from GMM method}
}
\description{
Estimate enveloped central space via global efficient estimator
}
\examples{
set.seed(500)
X=mvrnorm(500, mu = rep(0, 3), Sigma =diag(c(1,1,0.01)) , tol=1e-6)
Y=rbinom(500,1,1/(1+exp(-0.5*X[,1]-1*X[,2])))*(X[,1]+2*X[,2])
result=globaleff(X,Y,1)
Gammaest=result$Gamma1
Gammaest
}
