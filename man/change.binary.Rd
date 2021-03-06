% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChangePoint.R
\name{change.binary}
\alias{change.binary}
\title{Changepoint Model with Binary Discrete Distribution}
\usage{
change.binary(v, t = seq(1, length(v)), x = rep(NA, length(v)), pre = NA)
}
\arguments{
\item{v}{Vector that contains response variable of interest.}

\item{t}{Vector that contains threshold variable over which to study change
(ex. Time, outcome count in succession); default is sequence from 1 to
total number of observations.}

\item{x}{Vector that contains hidden parameter, should you wish to test an
interaction of original data with another variable (i.e., a hierarchical
model); NA by default.}

\item{pre}{A value represents preselected changepoint for any model, should
the user wish to test a specific point in data; NA by default.}
}
\value{
Table that contains location of changepoint (either preselected
  or determined rigorously) and its associated significance value.
}
\description{
Fit changepoint model on data to test if there is a significant change in
binary probability of one outcome versus other outcome.
}
\examples{
r = rnorm(400,0,1) * 10
r1 = r[1:200]
r2 = r[201:400]
a = rbinom(200,1,
           exp(r1)/(1+exp(r1)))
b = rbinom(200,1,
           exp(r2+1)/(1+exp(r2+2)))
d = data.frame(x=r,y=c(a,b),t=1:400)

change.binary(v = d$y, t= d$t, x = d$x)
change.binary(v = d$y, t= d$t, pre = 200)

}
