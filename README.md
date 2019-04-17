# Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB)

This Shiny App will allow you to explore the Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB), a database that integrates and compares the transcriptional, chromatin accessibility, and DNA methylation landscapes of vascular endothelial cells (ECs). Unbiased RNA-seq, ATAC-seq, and MethylC-seq enabled exploration of these landscapes among postnatal day 7 (P7) vascular ECs isolated from brain, liver, lung, and kidney of Tie2-GFP transgenic mice. These datasets should provide a foundation to determine the factors that are associated with EC heterogeneity. The app also allows exploration of EC gene expression in the developing CNS at single cell resolution.

# Running the app locally

To use the app on your local computer, simply run the following commands from an R console:

```R
install.packages("shiny")
library(shiny)
runGitHub("EC_Genomics","mfsabbagh")
```

# Citation

Please cite as: Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. "Transcriptional and epigenomic landscapes of CNS and non-CNS vascular endothelial cells." eLife 2018;7:e36187 DOI: 10.7554/eLife.36187
