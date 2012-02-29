\name{colon}
\alias{colon}
\docType{data}
\title{expression levels of colon tissues for human genes.}
\description{
 The \code{colon} includes data expression levels of 40 tumor and 22 normal colon tissues for 6,500 human genes. We here select of 2,000 genes with highest minimal intensity across the samples.
}
\usage{data(colon)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x}{a matrix with 62 rows and 2000 columns}
    \item{y}{class labels}
  }
}
\references{
Alon, U., Barkai, N., Notterman, D., Gish, K., Ybarra, S., Mack, D., and Levine, A. (1999), "Broad patterns of gene expression revealed by clustering analysis of tumor and normal colon tissues probed by oligonucleotide arrays,"
\emph{Proceedings of the National Academy of Sciences, 96, 6745.}
}
\examples{
data(colon)
}
\keyword{datasets}
