\name{penaltyParameter}
\alias{penaltyParameter}
\title{Compute the penalty parameter for the model.}
\description{Find the best penalty parameter \eqn{C} for the generalized distance weighted discrimination (DWD) model.}
\usage{penaltyParameter(X,y,expon,rmzeroFea = 1, scaleFea = 1)}
\arguments{
	\item{X}{A \eqn{d} x \eqn{n} matrix of \eqn{n} training samples with \eqn{d} features.}
	\item{y}{A vector of length \eqn{n} of training labels. The element of \code{y} is either -1 or 1.}
	\item{expon}{A positive number representing the exponent \eqn{q} of the residual \eqn{r_i} in the generalized DWD model. Common choices are \code{expon = 1,2,4}.}
	\item{rmzeroFea}{Switch for removing zero features in the data matrix. Default is set to be 1 (removing zero features).}
	\item{scaleFea}{Switch for scaling features in the data matrix. This is to make the features having roughly similar magnitude. Default is set to be 1 (scaling features).}
}
\details{
 The best parameter is empirically found to be inversely proportional to the typical distance between different samples raised to the power of (\eqn{expon+1}). 
 It is also dependent on the sample size \eqn{n} and feature dimension \eqn{d}.
}
\value{A number which represents the best penalty parameter for the generalized DWD model.}
\references{
Lam, X.Y., Marron, J.S., Sun, D.F., and Toh, K.C. (2018)
  ``Fast algorithms for large scale generalized distance weighted discrimination", \emph{Journal of Computational and Graphical Statistics}, forthcoming. \cr
  \url{https://arxiv.org/abs/1604.05473}
}
\author{Xin-Yee Lam, J.S. Marron, Defeng Sun, and Kim-Chuan Toh}
\keyword{DWD}

\examples{
# load the data
data("mushrooms")
# calculate the best penalty parameter
C = penaltyParameter(mushrooms$X,mushrooms$y,expon=1)
}
