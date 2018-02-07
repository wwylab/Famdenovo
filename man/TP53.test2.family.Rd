\name{TP53.test2.family}
\alias{TP53.test2.family}
\docType{data}
\title{
Family information data
}
\description{
Family information data used in function Famdenovo
}
\usage{TP53.test2.family}
\format{
  A data frame containing 5 variables.
  \describe{
    \item{\code{id}}{individual id}
    \item{\code{fid}}{father id}
    \item{\code{mid}}{mother id}
    \item{\code{gender}}{gender. 0: female, 1: male}
    \item{\code{age}}{If the individual was alive at the last contact date, use the age at that time. otherwise, use the age at death.}
  }
}
\examples{
data(TP53.test2.family)
TP53.test2.family
}
\keyword{datasets}
