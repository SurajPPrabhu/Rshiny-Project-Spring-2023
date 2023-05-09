# Rshiny-Project-Spring-2023 for the bF 591 class 
# mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals

This is a web application developed using R Shiny for exploring mRNA-Seq expression profiling data from the human post-mortem BA9 brain tissue for Huntington’s Disease (HD) and neurologically normal (NN) individuals. The data used in this app is from the paper "Transcriptomic profiling of post-mortem human brains reveals commonalities in gene expression changes in Huntington's disease and other neurodegenerative diseases" by Johnson et al. (2016).

The app allows users to visualize and explore the differential expression of genes between HD and NN samples, as well as the expression patterns of individual genes in each sample. The app includes several tabs, each with different functionality:

- **Data:** Displays a table of the original data used in the app. The table can be filtered and sorted by column.
- **Counts:** Displays a table of gene expression counts for each sample, as well as summary statistics and a dropdown menu for selecting a gene of interest.
- **PCA:** Displays a principal component analysis (PCA) plot of the gene expression data, with the ability to color points by HD or NN status.
- **DE Analysis:** Performs differential expression analysis between HD and NN samples, and displays a table of the results. Users can adjust the p-value and log-fold change thresholds for the analysis.
- **Volcano Plot:** Displays a volcano plot of the differential expression results, with the ability to adjust the p-value and log-fold change thresholds for highlighting significant genes.
- **Heatmap:** Displays a heatmap of the gene expression data, with the ability to cluster genes and samples by similarity.
- **Density Plots:** Displays density plots of the gene expression data, with one plot for each numeric column in the data.

To use the app, users can download the `app.R` file and run it in RStudio, or access the app online at [insert link to shinyapps.io].

## Acknowledgements

This app was developed as a final project for BF591 at Boston University, taught by Dr. Paige Nong.
