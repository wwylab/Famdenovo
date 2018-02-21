# Famdenovo
Calculate the probability of the de novo status for a germline mutation in familial diseases.

The inputs are family data, cancer data and mutation data and the person(s) you want to analyze. All individuals should from one family. The output is the probability of any TP53 mutation being de novo, one TP53 mutation carrier per line. Each line contains three elements: "family id", "individual id" and "predicted denovo probability", respectively.

# Installation
To install the package, start R and enter:

devtools::install_github("wwylab/Famdenovo/Famdenovo_0.1.1")

For more information, please visit:

https://bioinformatics.mdanderson.org/main/Famdenovo
