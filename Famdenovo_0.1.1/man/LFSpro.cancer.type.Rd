\name{LFSpro.cancer.type}
\alias{LFSpro.cancer.type}
\docType{data}
\title{
Predefined cancer types and the corresponding number in LFSpro
}
\description{
We classified the cancers into 11 groups according to NCCCN Guidelines Version 1.2012 Li-Fraumeni Syndrome criteria. And then we coded these different groups of cancer into different number. 

sts(soft tissue sarcoma): 1

ost(osterosarcoma): 2

brain(brain tumor): 3

breast(breast cancer): 4

acc(adrenocortical carcinoma): 5

leukemia: 6

lung(lung bronchoalveolar cancer): 7

choroid(choroid plexus carcinoma): 8

other.lfs(other LFS spectrum cancers): 20

non.lfs(non LFS spectrum invasive cancers): 50

benign(benign tumors): 100
}
\usage{LFSpro.cancer.type}
\format{
A vector with names of different cancer types.

        sts       ost     brain    breast       acc  leukemia      lung   choroid other.lfs   non.lfs benign  
        
        1         2         3         4         5         6         7         8        20        50
      100 
}

\examples{
data(LFSpro.cancer.type)
LFSpro.cancer.type
}
\keyword{datasets}
