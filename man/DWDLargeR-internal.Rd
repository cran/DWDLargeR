\name{DWDLargeR-internal}
\title{Internal DWDLargeR functions}
\alias{sigma_update}
\alias{linsysolve}
\alias{smw}
\alias{fnorm}
\alias{polyRootsNewton}
\alias{psqmr}
\alias{precondfun}
\alias{vecMultiply}
\description{Internal DWDLargeR functions.}
\usage{
sigma_update(iter)
linsysolve(R,r)
smw(R,Z,ZT,yunit,schurmat,schurvec,normy,const,r)
fnorm(x)
polyRootsNewton(c,q,sigma,x0)
psqmr(ff,b,L,x0,tol,maxit)
precondfun(L,r)
vecMultiply(ff,x)
}
\author{Xin-Yee Lam, J.S. Marron, Defeng Sun, and Kim-Chuan Toh}
\references{
Lam, X.Y., Marron, J.S., Sun, D.F., and Toh, K.C. (2018)
  ``Fast algorithms for large scale generalized distance weighted discrimination", \emph{Journal of Computational and Graphical Statistics}, forthcoming. \cr
  \url{https://arxiv.org/abs/1604.05473}
}

\details{These are not intended for use by users.}
\keyword{internal}
