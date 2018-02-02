\name{TP53.test1.mutation}
\alias{TP53.test1.mutation}
\docType{data}
\title{
Mutation information data
}
\description{
Mutation information data used in function Famdenovo
}
\usage{TP53.test1.mutation}
\format{
  A data frame containing 2 variables.
  \describe{
    \item{\code{id}}{individual id}
    \item{\code{mut.state}}{mutation status of the person. "W" - wild type; "M" - mutated. Individuals who are not sequenced are not included in the mutation information data.}
  }
}
\examples{
data(TP53.test1.mutation)
TP53.test1.mutation
}
\keyword{datasets}
