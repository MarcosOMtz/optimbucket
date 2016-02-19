\name{auc.wroc}
\alias{auc.wroc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
AUC and Gini Calculation
}
\description{
Calculate predictive power of a variable using a wroc object.
}
\usage{
auc.wroc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class 'wroc'.}
}
\details{
The ROC curve here in fact uses True Negative Rate (tnr) vs False Negative Rate (fnr) instead of True Positive Rate (tpr) vs False Positive Rate (fpr) to make it easier to work in credit scoring. The measures are equivalent.
}
\value{
An object of class 'wauc' including several variations on the same performance measure:
\item{auc}{Area under the ROC curve.}
\item{gini_two_way}{abs(2*auc - 1)}
\item{gini_up}{Twice the area above y=x and the ROC curve.}
\item{gini_down}{Twice the area below y=x and the ROC curve.}
\item{gini}{The larger between gini_up and gini_down.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Felipe Gerard
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
See ?wroc for examples.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~AUC }
\keyword{ ~Gini }
\keyword{ ~ROC }
