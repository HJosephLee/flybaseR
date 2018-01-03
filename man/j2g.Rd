\name{j2g}
\alias{j2g}
\title{Open FlyBase gene reports using a web browser}
\usage{
j2g(x, n = 10)
}
\arguments{
\item{x}{a character or vector. FlyBase IDs or symbols.}

\item{n}{The maximum number of genes to be opened.  Default = 10.}
}
\description{
The function takes FlyBase IDs (e.g. FBgn#######) or Gene symbols as an input, and opens FlyBase gene report from the default web browser (e.g. Jump2Gene).
}
\examples{
j2g("e2f")
j2g(c("ovo", "otu"))
j2g("FBgn0000504")
j2g(c("FBgn0000504", "tra"))
}
\keyword{flybase}
