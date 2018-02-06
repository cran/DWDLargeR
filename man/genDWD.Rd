\name{genDWD}
\alias{genDWD}
\title{Solve the generalized distance weighted discrimination (DWD) model.}
\description{Solve the generalized DWD model by using a symmetric Gauss-Seidel based alternating direction method 
  of multipliers (ADMM) method.}
\usage{genDWD(X,y,C,expon, tol = 1e-5, maxIter = 2000, method = 1, printDetails = 0,
             rmzeroFea = 1, scaleFea = 1)}
\arguments{
	\item{X}{A \eqn{d} x \eqn{n} matrix of \eqn{n} training samples with \eqn{d} features.}
	\item{y}{A vector of length \eqn{n} of training labels. The element of \code{y} is either -1 or 1.}
	\item{C}{A number representing the penalty parameter for the generalized DWD model.}
	\item{expon}{A positive number representing the exponent \eqn{q} of the residual \eqn{r_i} in the generalized DWD model. Common choices are \code{expon = 1,2,4}.}
	\item{tol}{The stopping tolerance for the algorithm. (Default = 1e-5)}
	\item{maxIter}{Maximum iteration allowed for the algorithm. (Default = 2000)}
	\item{method}{Method for solving generalized DWD model. The default is set to be 1 for the highly efficient sGS-ADMM algorithm. User can also select \code{method = 2} for the directly extended ADMM solver.}
	\item{printDetails}{Switch for printing details of the algorithm. Default is set to be 0 (not printing).}
	\item{rmzeroFea}{Switch for removing zero features in the data matrix. Default is set to be 1 (removing zero features).}
	\item{scaleFea}{Switch for scaling features in the data matrix. This is to make the features having roughly similar magnitude. Default is set to be 1 (scaling features).}
}
\details{
 This is a symmetric Gauss-Seidel based alternating method of multipliers (sGS-ADMM) algorithm for solving the generalized DWD model of the following formulation:
 \deqn{\min \sum_i \theta_q (r_i) + C e^T x_i} subject to the constraints \deqn{Z^T w + \beta y + \xi - r = 0, ||w||<=1, \xi>=0,} \cr
 where \eqn{Z = X diag(y)}, \eqn{e} is a given positive vector such that \eqn{||e||_\infty = 1}, and \eqn{\theta_q}
 is a function defined by \eqn{\theta_q(t) = 1/t^q} if \eqn{t>0} and \eqn{\theta_q(t)=\infty} if \eqn{t<=0}.
}

\value{A list consists of the result from the algorithm.
	\item{w}{The unit normal of hyperplane that distinguishes the two classes.}
	\item{beta}{The distance of the hyperplane to the origin (\eqn{\beta} in the above formulation).}
	\item{xi}{A slack variable of length \eqn{n} for the possibility that the two classes may not be separated cleanly by the hyperplane (\eqn{\xi} in the above formulation).}
	\item{r}{The residual \eqn{r:= Z^T w + \beta y + \xi}.}
	\item{alpha}{Dual variable of the generalized DWD model.}
	\item{info}{A list consists of the information from the algorithm.}
	\item{runhist}{A list consists of the run history throughout the iterations.}
}

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
# solve the generalized DWD model
result = genDWD(mushrooms$X,mushrooms$y,C,expon=1)
}
