\name{id.converter}
\alias{id.converter}
\title{Converting or Updating FlyBase IDs or Symbols}
\usage{
id.converter(x, symbols = F, bundle.size = 1000, DmelOnly = T, polite.access = 0, diehard.symbols = F)
}
\arguments{
\item{x}{a vector. FlyBase IDs or names to be converted.}

\item{symbols}{Logical.  If TRUE, the output will be gene symbols, rather than FlyBase IDs.  Default = F}

\item{bundle.size}{Numeric.  The number of FlyBase IDs or symbols to be submitted to FlyBase at once. Default is 1,000 if there are less than 100 symbols; 100 if more than 1,000 symbols.  Reduce the number down if Timeout error occurs.}

\item{DmelOnly}{Logical.  If TRUE, non-melanogaster gene IDs will be ignored.  Default = T.}

\item{polite.access}{Numeric.  Intervals between FlyBase access for each bundle as seconds.  Default = 0.}

\item{diehard.symbols}{Logical.  If TRUE, ntervals between FlyBase access for each bundle as seconds.  Default = 0.}

\item{convert.into}{"genes", "transcripts", or "polypeptides". "g", "t", or "p" is also possible. If missing, the IDs will be updated to the most recent IDs only.}
}
\description{
The function takes FlyBase IDs (e.g. FBgn0000003) or Gene symbols as an input, and converts it into updated IDs or symbols using the FlyBase ID converter (web).
The function accesses FlyBase, so requires internet-connection.  FlyBase ID inputs are bundled as 1,000.  100 for symbols.
FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::).  Genes that does not have matching IDs will be shown as "unknown".
Certain gene symbols would appear as "unknown" even if the gene exists, and have FlyBase IDs. This is because the ID converter in FlyBase website cannot convert the gene.  For example, CG31976 cannot be converted by FlyBase, although you can find the gene from the gene report.  Setting diehard.symbols = T will look for gene report pages of such unconvertible genes one by one.  The process is essentially slow because it accesses FlyBase for each gene.
}
\examples{
id.converter(x, symbols = T)
id.converter(x, bundle.size = 50, be.polite = 10, convert.into = "transcripts")
id.converter(x, symbols = T, bundle.size = 50, diehard.symbols = T)
}
\keyword{flybase}
  
