\name{wroc}
\alias{wroc}
\alias{wroc.default}
\alias{optimize.wroc}
\alias{subset.wroc}
\alias{plot.wroc}
\alias{print.wroc}

\title{
Create and handle wroc objects
}
\description{
Functions for creating and handling wroc objects. You can get optimal bucketings with optimal.wroc and join buckets manually with subset.wroc.
}
\usage{
function(predictions, labels, ngroups=50, level.bad=1, col.bad=1)

\method{wroc}{default}(predictions, labels, ngroups=50, level.bad=1, col.bad=1)
\method{print}{wroc}(x, ...)
\method{plot}{wroc}(x, type = c('accum','roc','trend','woe'))
\method{optimize.wroc}{wroc}(x, trend = c('auto','upper','lower'))
\method{subset}{wroc}(x, buckets = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictions}{Either a numeric vector or an ordered factor.}
  \item{labels}{Either a vector with values for the predictor variable or a 2-column matrix containing counts for each class for each value of predictions.}
  \item{ngroups}{Precision to use for the ROC curve. Similar to rounding idea.}
  \item{level.bad}{If labels is a vector, the level corresponding to the positive class.}
  \item{col.bad}{If labels is a matrix, the column containing the counts for the positive class.}
  \item{x}{An object of type wroc.}
  \item{type}{The type of plot to show. 'accum' is similar to a ROC curve, but using cumulative positive and negative examples proportions, where 'roc' uses traditional right-handed proportions. 'trend' shows the trend of the probability of the positive class. 'woe' whos the trend of the Weight of Evidence (WoE), given by 'log((G_i/sum(G_i))/(B_i/sum(B_i)))'}
  \item{trend}{Trend to use for calculating Gini index. Can be one of 'upper', 'lower' or 'auto'. 'auto' chooses the one among the other two with the most points below the y=x line on ROC space.}
  \item{buckets}{Bucket numbers to join to the adjacent right bucket.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class `wroc`, `optimal.wroc` or `subset.wroc`, with the following elements:
\item{info}{A data.frame with most of the information about the bucketing, including counts of positive / negative classes, probabilities and limits of the original variable.}
\item{ngroups}{The actual number of groups the bucketing has, equal to the number of rows of `info`, excluding the bucket number 0.}
\item{call}{The original call to `wroc`.}
\item{jumps}{(optimize.wroc only, for internal use.) The amount of buckets to the right the optimization skipped. Read from left to right.}
\item{groups}{New buckets used for summarization.}
\item{call.optimize}{(optimize.wroc only) The call to optimize.wroc}
\item{call.subset}{(subset.wroc only) The call to subset.wroc}
}
\references{
optimize.wroc is based on the original idea by Ramón Jiménez.
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

## Load required package
library(ggplot2)
library(optimbucket)

## Generate data
ng <- 10000
nb <- 1000

goods <- data.frame(
  x = rnorm(ng, mean = ifelse(runif(ng) > 0.2, 10, -10)),
  y = 0
)
bads <- data.frame(
  x = rnorm(nb),
  y = 1
)

d <- rbind(goods, bads) %>%
  mutate(y = factor(y))

## Plot positive and negative distributions
ggplot(d, aes(x, fill=y)) +
  geom_density(alpha=0.5)

## Simple ROC object
wr <- wroc(d$x, d$y, ngroups=30)
auc.wroc(wr)
plot(wr)
plot(wr, type='trend')
plot(wr, type='woe')

## Optimal ROC object with automatic trend
wr2 <- optimize.wroc(wr, 'auto')
auc.wroc(wr2)
plot(wr2)
plot(wr2, type='trend')
plot(wr2, type='woe')

## Manual bucketing
wr3 <- subset(wr, c(1:7,10:29))
auc.wroc(wr3)
plot(wr3)
plot(wr3, type='woe')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~ROC }
\keyword{ ~bucket }
\keyword{ ~credit }
\keyword{ ~scoring }
