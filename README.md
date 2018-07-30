# EC_Genomics

This Shiny App will allow you to explore the Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB), a database that integrates and compares the transcriptional, chromatin accessibility, and DNA methylation landscapes of vascular endothelial cells. Unbiased RNA-seq, ATAC-seq, and MethylC-seq enabled exploration of these landscapes among four populations of early postnatal murine ECs to determine the factors that are associated with EC heterogeneity. The app also allows exploration of endothelial cell gene expression in the developing CNS at single cell resolution.

The app is also accessible via https://markfsabbagh.shinyapps.io/vectrdb/

# Running the app locally

To use the app on your local computer, simply run the following commands from an R console:

```R
install.packages("shiny")
library(shiny)
runGitHub("EC_Genomics","mfsabbagh")
```

# Citation

Please cite: Mark F. Sabbagh, Jacob S. Heng, Chongyuan Luo, Rosa G. Castanon, Joseph R. Nery, Amir Rattner, Loyal A. Goff, Joseph R. Ecker, and Jeremy Nathans. Vascular Endothelial Cell Trans-omics Resource Database (VECTRDB). https://markfsabbagh.shinyapps.io/vectrdb/
