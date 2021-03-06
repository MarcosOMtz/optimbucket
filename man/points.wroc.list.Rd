% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R, R/wroc.R
\name{point_transform}
\alias{point_transform}
\alias{points.wroc.list}
\title{Transformation to Points}
\usage{
point_transform(point.transform, woe, var = NULL)

\method{points}{wroc.list}(x, model, base.points = 600, base.odds = 1/20,
  pdo = 15, point.decimals = 0, reoptimize = TRUE,
  reoptimize.point.tol = 0, prefix.to.remove = "woe_")
}
\arguments{
\item{x}{An object of class \code{wroc.list}}

\item{model}{An object of class \code{glm}. Note that the variables used must
be included in \code{x}}

\item{base.points}{Points to be given at \code{base.odds}}

\item{pdo}{Points to double the odds}

\item{point.decimals}{Precision parameter for the points}

\item{reoptimize}{Should buckets with less than \code{reoptimize.point.tol}
steps be pasted together?}

\item{reoptimize.point.tol}{See \code{reoptimize}. If \code{reoptimize =
FALSE}, this does nothing}

\item{prefix.to.remove}{The prefix to remove from the \code{glm}'s variable
names (since usually the \code{wroc} objects have var1, var2, ... but the
glm will have woe_var1, woe_var2, ...)}
}
\value{
An object of class \code{points.wroc} (apart from its other classes)
  which contains a \code{points} variable in its \code{info} component and a
  \code{point.transform} component containing all the relevant parameters.
}
\description{
Given a wroc.list and a logistic regression object, calculates the points for
each level of each variable and yields a list of points.wroc object.
}

