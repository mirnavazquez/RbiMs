

# rbims <img src="man/figures/Logo-rRbiMs.png" alt="Logo Rbims" width="150" align="right"/>

[![R-CMD-check](https://github.com/mirnavazquez/RbiMs/workflows/R-CMD-check/badge.svg)](https://github.com/mirnavazquez/RbiMs/actions)

**rbims** (Reconstruction of Bin Metabolisms) is an R package designed to streamline the functional analysis and visualization of metagenome-assembled genomes (MAGs). It supports annotation integration from KEGG, InterProScan, dbCAN, and MEROPS, allowing researchers to quantify gene presence, abundance, and pathway coverage across microbial genomes.

The package includes a curated database for hydrocarbon degradation pathways (aerobic and anaerobic) and provides tools to generate publication-ready visualizations such as heatmaps and bubble plots. It is designed to assist in exploratory trait analysis and early-stage hypothesis generation in genome-resolved metagenomics.

---

## ✨ Features

- Import functional annotations from KEGG, InterProScan, dbCAN, MEROPS, and PICRUSt2
- Calculate presence, abundance, and pathway coverage per MAG
- Subset data by gene, enzyme, pathway, or domain
- Visualize functional traits with customizable bubble plots and heatmaps
- Integrate sample or genome metadata
- Export data frames and visualizations for publication
- Curated database for hydrocarbon degradation not covered in KEGG

---

## 🚀 Quick Install

```r
install.packages("devtools")
library(devtools)
install_github("mirnavazquez/RbiMs")
library(rbims)
```

If you are on **macOS**, install [XQuartz](https://www.xquartz.org/).  
If using **Ubuntu**, install system dependency: `libcairo2-dev`.

---

## 🧬 Case Study: Oil-Enriched Marine MAGs

A complete example using MAGs from a hydrocarbon enrichment experiment is available in the folder [`/Hidrocarburos`](https://github.com/mirnavazquez/RbiMs/tree/main/Hidrocarburos), including annotation files and code to reproduce the figures in our manuscript.


---

## 👩‍💻 Contributors

- **Mirna Vázquez-Rosas-Landa** – lead developer  
- **Karla P. López-Martínez** – co-developer, documentation, manuscript  
- **Stephanie Hereira-Pacheco** – functions and documentation  
- **Diana Hernández-Oaxaca** – conda environment setup  
- **Frida López-Ruiz** – testing and documentation  

---

## 📚 References

- Kanehisa M. and Goto S., Nucleic Acids Res. (2000)  
- Kanehisa M. et al., Protein Sci. (2019)  
- Kanehisa M. et al., Nucleic Acids Res. (2021)  
- [DiTing – hydrocarbon cycles definitions](https://github.com/xuechunxu/DiTing)

---

## 🌐 Website

Full documentation: [https://mirnavazquez.github.io/RbiMs](https://mirnavazquez.github.io/RbiMs)

