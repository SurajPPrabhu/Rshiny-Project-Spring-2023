# Rshiny-Project-Spring-2023 for the bF 591 class 
# RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression

This is a web application developed using R Shiny for exploring of the data from the paper titled:"RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression" by Labadorf et al.

The app allows users to visualize and explore the differential expression of genes between HD and NN samples, as well as the expression patterns of individual genes in each sample. The app includes several tabs, each with different functionality:

- **Data:** Displays a table of the original data used in the app. The table can be filtered and sorted by column.
- **Counts:** Displays a table of gene expression counts for each sample, as well as summary statistics and a dropdown menu for selecting a gene of interest.
- **PCA:** Displays a principal component analysis (PCA) plot of the gene expression data, with the ability to chose which Principal components the user wants.
- **DE Analysis:** Performs differential expression analysis between HD and NN samples, and displays a table of the results. Users can adjust the p-value and log-fold change thresholds for the analysis.
- **Volcano Plot:** Displays a volcano plot of the differential expression results, with the ability to adjust the p-value and log-fold change thresholds for highlighting significant genes.
- **Heatmap:** Displays a heatmap of the gene expression data, with the ability to cluster genes and samples by similarity.
- **Density Plots:** Displays density plots of the gene expression data, with one plot for each numeric column in the data.

To use the app, users can download the `app.R` file and run it in RStudio.

## Acknowledgements

This app was developed as a final project for BF591 at Boston University, taught by Dr. Adam

## Contents
- app- contains the app.r.
- Data- Contains the data used to test and build this app.
