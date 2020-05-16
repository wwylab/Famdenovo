# Famdenovo
This is a software that predicts the de novo status of deleterious germline mutation(s) in Mendelian diseases based on family history.  

The input data include family pedigree, cancer type, mutation information, and the person(s) for whom the probabilities will be computed.  The output is the probability of the deleterious mutation(s) being de novo. Famdenovo has been validated for predicting the de novo status of deleterious TP53 mutations. The validation study is described in Fan et al. 2020 at https://www.biorxiv.org/content/10.1101/2020.02.10.942409v1. 

# Installation
To install the package, start R and enter:

devtools::install_github("wwylab/Famdenovo/Famdenovo_0.1.1")

For more information, please visit:

https://bioinformatics.mdanderson.org/main/Famdenovo
