SISPA: Method for Sample Integrated Set Profile Analysis

Authors:
    - Bhakti Dwivedi (1)
    - Jeanne Kowalski (1,2)
  
Maintainer:
    - Bhakti Dwivedi <bhakti.dwivedi@emory.com>

Affiliations:
    (1) Winship Cancer Institute, Emory University, Atlanta, GA, 30322, USA
    (2) Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University, Atlanta, GA 30322, USA

Installation:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version='devel')
BiocManager::install("SISPA")


Getting Started:
    See ./vignettes/SISPA.pdf or ./vignettes/SISPA.html

To cite the SISPA method: 
Kowalski J, Dwivedi B, Newman S, Switchenko JM, Pauly R, Gutman DA, Arora J, Gandhi K, Ainslie K, Doho G, Qin Z, Moreno CS, Rossi MR, Vertino PM, Lonial S, Bernal-Mizrachi L, Boise LH. Gene integrated set profile analysis: a context-based approach for inferring biological endpoints. Nucleic Acids Res. 2016 Apr 20;44(7):e69. doi: 10.1093/nar/gkv1503. Epub 2016 Jan 29. PubMed PMID: 26826710; PubMed Central PMCID: PMC4838358

To cite the SISPA Bioconductor package:
Dwivedi B and Kowalski J (2016). SISPA: SISPA: Method for Sample Integrated Set Profile Analysis. R package version 1.8.0. 

To cite the SISPA shiny web-tool (http://shinygispa.winship.emory.edu/shinySISPA/):
Dwivedi B and Kowalski J. shinySISPA: A web tool for defining sample groups using gene sets from multiple-omics data. F1000Research 2018, 7:213 (doi: 10.12688/f1000research.13934.1) 