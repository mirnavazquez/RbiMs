 <!-- badges: start -->
  [![R-CMD-check](https://github.com/mirnavazquez/RbiMs/workflows/R-CMD-check/badge.svg)](https://github.com/mirnavazquez/RbiMs/actions)
  <!-- badges: end -->
  
# Rbims

R tools for reconstructing bin metabolisms.

## Quick install

In R terminal:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("KEGGREST", "devtools"))
```

```
library(devtools)
install_github("mirnavazquez/RbiMs")
library(rbims)
```
  
## Overview 

![](inst/rRbiMs-3.png)

* Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000).
* Kanehisa, M; Toward understanding the origin and evolution of cellular organisms. Protein Sci. 28, 1947-1951 (2019).
* Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M., and Tanabe, M.; KEGG: integrating viruses and cellular organisms. Nucleic Acids Res. 49, D545-D551 (2021).
* Rawlings, N.D., Waller, M., Barrett, A.J. & Bateman, A. (2014) MEROPS: the database of proteolytic enzymes, their substrates and inhibitors. Nucleic Acids Res 42, D503-D509.


