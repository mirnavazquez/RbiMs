
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **Rbims** <img src="man/figures/Logo-rRbiMs.png"  width="150" height="150" align="right" />

<!-- badges: start -->
<!-- badges: end -->

R tools for reconstructing bin metabolisms.

## Quick install

In R terminal:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("KEGGREST", "devtools"))
```

And the development version from
[GitHub](https://github.com/mirnavazquez/RbiMs) with:

``` r
library(devtools)
install_github("mirnavazquez/RbiMs")
library(rbims)
```

## Overview

![](inst/rRbiMs-3.png)

-   Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and
    Genomes. Nucleic Acids Res. 28, 27-30 (2000).
-   Kanehisa, M; Toward understanding the origin and evolution of
    cellular organisms. Protein Sci. 28, 1947-1951 (2019).
-   Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M., and
    Tanabe, M.; KEGG: integrating viruses and cellular organisms.
    Nucleic Acids Res. 49, D545-D551 (2021).
-   [DiTing](https://github.com/xuechunxu/DiTing) cycles definition.
