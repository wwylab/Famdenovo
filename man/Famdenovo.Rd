\name{Famdenovo}
\alias{Famdenovo}
\title{
  Calculate Pr(denovo) for a TP53 mutation
}
\description{
  Calculate the probability of a TP53 mutation being denovo.
}
\usage{
Famdenovo(family, cancer, mutation, person.id, gene = TP53)
}
\arguments{

\item{family}{The input should be a data frame. The family data should include the following columns with the corresponding column names:

  id: index of the person. All individuals should from one family.

  fid: index of the person's father. If the individual is the founder of the pedigree, set it as NA.

  mid: index of the person's mother. If the individual is the founder of the pedigree, set it as NA.

  gender: gender of the person. 0 - female; 1 - male

  age: age of the person. If the individual is alive, it is as the current age. Otherwise, set it as the age of death.
}

\item{cancer}{The input should be a data frame. The cancer data should include the following columns with the corresponding column names:

  id: index of the person

  cancer.type: type of the cancer. We divided all the cancers into 11 groups according to NCCCN Guidelines Version 1.2012 Li-Fraumeni Syndrome criteria. Check "LFSpro.cancer.type" for details.

  diag.age: The age when the individual was diagnosed with cancer.
}

\item{mutation}{The input should be a data frame. The mutation data should include the following columns with the corresponding column names:

  id: index of the person

  mut.state: mutation status of the person. "W" - wild type; "M" - mutated. Individuals who are not sequenced are not included in the mutation information data. 
}

\item{person.id}{The input should be either character string(s) or numrical value(s) of the person(s) you want to analyze.}

\item{gene}{The input should be character string(s). The default value is "TP53". We will add other genes in the future. .}
}

\details{
The inputs are family data, cancer data and mutation data and the person(s) you want to analyze. All individuals should from one family. The output is the probability of any TP53 mutation being de novo, one TP53 mutation carrier per line. Each line contains three elements: "family id", "individual id" and "predicted denovo probability", respectively. 
}

\references{
Peng G, Bojadzieva J, Ballinger ML, et al: Estimating TP53 Mutation Carrier Probability in Families with Li-Fraumeni Syndrome Using LFSPRO. Cancer Epidemiology, Biomarkers & Prevention, 2017

Xuedong Pan, Elissa B. Dodd, Fan Gao, et al: A general probabilistic algorithm to predict de novo mutations in familial diseases as demonstrated in Li-Fraumeni Syndrome, Statistics in Oncology, 2018
}

\examples{
library(Famdenovo)

# Test case 1
data(TP53.test1.family)
data(TP53.test1.cancer)
data(TP53.test1.mutation)
person.id <- c(1, 3, 1001, 1002, 1003)
output1 <- Famdenovo(TP53.test1.family, TP53.test1.cancer, person.id, TP53.test1.mutation,); output1

# Test case 2
data(TP53.test2.family)
data(TP53.test2.cancer)
person.id <- c("SMN-1186-0000", "SMN-1186-0702", "SMN-1186-0703", "SMN-1186-0701")
output2 <- Famdenovo(TP53.test2.family, TP53.test2.cancer, person.id); output2

# Compare with known denovo/familial status
data(TP53.test2.counselees)
valid.person.id <- person.id[person.id \%in\% unique(TP53.test2.counselees$id)]
TP53.test2.counselees[TP53.test2.counselees$id == valid.person.id, c("id", "cou.state")]
merge(TP53.test2.counselees[TP53.test2.counselees$id == valid.person.id, c("id", "cou.state")], output2, all.x=T)
}
