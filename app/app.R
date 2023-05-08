library(shiny)
library(shinydashboard)
library(tidyverse)
library(scales)
library(gplots)
library(RColorBrewer)
#defining the layout and sidebars
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Sample Information Exploration", tabName = "data"),
    menuItem("Counts", tabName = "counts")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "data",
            h2("Explore your data and its features!"),
            fluidRow(
              column(8, fileInput("file1", "Upload CSV file"))
            ),
            uiOutput("data_tabs")
    ),
    
    tabItem(tabName = "counts",
            h2("Counts data analysis"),  
            wellPanel(
              fluidRow(
                column(width = 7,
                       box(
                         title = "Upload CSV file containing counts information",
                         fileInput("file2", ""),
                         br(),
                         sliderInput("gvar", "Select the minimum percentile of variance desired in the genes",
                                     min = -50, max = 100, value = 50),
                         sliderInput("nz", "Select the minimum number of nonzero samples",
                                     min = 0, max = 100, value = 50),
                         actionButton("submit_btn", "Submit")
                       )
                )
              ),
              uiOutput("counts_tabs")
            )
    )
  )
)


#defining server side functions
server <- function(input, output) {
  options(shiny.maxRequestSize = 30*1024^2)
  
  # Read data file
  data_summary <- reactive({
    if (is.null(input$file1)) {
      return(NULL)
    } else {
      read.csv(input$file1$datapath, header = TRUE)
    }
  })
  
  # Define data_counts reactive function to read counts file
  data_counts <- reactive({
    if (is.null(input$file2)) {
      return(NULL)
    } else {
      read_csv(input$file2$datapath, col_names = TRUE,)
    }
  })
  
  # Define summary_data function to process counts data
  summary_data <- function(data_counts, gvar, nz_genes) {
    # Filter genes based on variance and number of non-zero samples
    genes_var_filtered <- data_counts[rowSums(data_counts > 0) >= nz_genes & apply(data_counts, 1, var) >= gvar, ]
    
    # Get summary statistics for the filtered data
    num_rows_filtered <- nrow(genes_var_filtered)
    num_rows_not_filtered <- nrow(data_counts) - num_rows_filtered
    pgene <- percent((num_rows_filtered / nrow(data_counts)))
    fgene <- percent((num_rows_not_filtered) / nrow(data_counts))
    
    # Create summary data frame
    summary_df <- data.frame(Description = c("Number of Genes", "Number of Samples", "Number of Genes (Filtered):", "% of genes filtered", "Number of Genes not filtered:", "%genes not filtered:"), Value = c(nrow(data_counts), ncol(data_counts), num_rows_filtered, pgene, num_rows_not_filtered, fgene))
    
    return(summary_df)
  }
  
  
  
  # Render the data tabs UI
  output$data_tabs <- renderUI({
    if (is.null(input$file1)) {
      return(NULL)
    } else {
      tabsetPanel(
        tabPanel("Counts"),
        tabPanel("Sample Information"),
        tabPanel("Plot my data")
      )
    }
  })
  
  # Render the counts tabs UI
  output$counts_tabs <- renderUI({
    if (is.null(data_counts())) {
      return(NULL)
    } else {
      tabsetPanel(
        tabPanel("Counts Summary", dataTableOutput("summary_data")),
        tabPanel("Counts Diagnostic Plots", plotOutput("median_plot"), plotOutput("numzero")),
        tabPanel("Heatmap", plotOutput("heatmap_plot")),
        tabPanel("PCA" ,selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                 selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2")
                , plotOutput("plot_pca"))
        )
    }
  })
  
  # Scatter plot of median count vs variance
  output$median_plot<- renderPlot({
    if (!is.null(data_counts())) {
      plot_tib <- data_counts()
        plot_tib<-plot_tib%>%mutate(Median = apply(plot_tib[-1], MARGIN = 1, FUN = median), 
               Variance = apply(plot_tib[-1], MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$Variance, probs = input$gvar/100)   #calculate percentile
      plot_tib <- plot_tib %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE")) #sort genes by percentile
      #plot scatter plot
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, Variance))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        labs(title = 'Plot of Median vs Variance.', subtitle = "Genes filtered out are in red. X and Y axes are log-scaled.")+
        scale_y_log10()+
        scale_x_log10()+
        theme_bw()+
        theme(legend.position = 'bottom')
      return(scatter)}})

  output$numzero <- renderPlot({
    if (!is.null(data_counts())) {
      data <- data_counts()
      tot_samp <- ncol(data) - 1  #store original number of samples
      #make a plot tibble
      data <- data %>%   
        mutate(Median = apply(data[-1], MARGIN = 1, FUN = median)) %>% 
        na_if(0)
      data$no_zeros <- rowSums(is.na(data))  #make new col, with counts.
      data$filtered <- ifelse(data$no_zeros <= input$nz, "TRUE", "FALSE")
      #plot scatter plot
      cols <- c("FALSE" = "red", "TRUE" = "black")
      scatter <- ggplot(data, aes(Median, no_zeros)) +
        geom_point(aes(color = filtered), alpha = 0.75) +
        scale_color_manual(values = cols) +
        scale_x_log10() +
        labs(
          title = 'Plot of Median vs Number of Non-Zero genes',
          subtitle = "Genes filtered out are in red. X-axis is log scaled."
        ) +
        theme_bw() +
        ylab('Number of samples with zero count') +
        theme(legend.position = 'bottom')
      return(scatter)
    } else {
      return(NULL)
    }
  })
  

  # Define the output plot
  output$heatmap_plot <- renderPlot({
    
    # Get the count data from the data_counts() function
    counts_df <- data_counts()
    
    # Set all zero values to NA
    counts_df <- na_if(counts_df, 0)
    
    # Create a new column with the number of zeros for each row
    counts_df$zeros_count <- rowSums(is.na(counts_df))
    
    # Filter out rows with more zeros than the input value nz
    counts_df <- filter(counts_df, zeros_count <= input$nz)
    
    # Log10 transform the count values, excluding the gene names and zeros count columns
    counts_df <- log10(counts_df[,!colnames(counts_df) %in% c("gene", "zeros_count")])
    
    # Create a new column with the variance for each row
    plot_df <- counts_df %>% 
      mutate(variance = apply(counts_df, MARGIN = 1, FUN = var))
    
    # Calculate the percentile value based on the input gvar
    perc_val <- quantile(plot_df$variance, probs = input$gvar/100, na.rm = TRUE)
    
    # Filter out rows with variance less than the percentile value
    plot_df <- filter(plot_df, variance >= perc_val)
    
    # Generate the heatmap plot using the heatmap.2 function from the gplots package
    hmap <- heatmap.2(as.matrix(plot_df[-ncol(plot_df)]), scale = "row", col = brewer.pal(9, "YlOrRd"))
    
    # Return the plot
    return(hmap)
  })
  
  output$plot_pca <- renderPlot({
    if (is.null(data_counts())) {
      return(NULL)
    } else {
      # calculate PCA
      filt_tib <- data_counts() %>% 
        mutate(variance = apply(data_counts()[-1], MARGIN = 1, FUN = var), .after = gene) # calculate variance for filtering
      perc_val <- quantile(filt_tib$variance, probs = input$gvar/100, na.rm = TRUE)   # calculate percentile
      filt_tib <- filter(filt_tib, variance >= perc_val) # filter the tibble
      pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) # transpose the data and perform PCA
      
      # extract variance explained by each component
      var_explained <- summary(pca_res)$importance[2,] * 100
      
      # get selected principal components
      comp1 <- which(colnames(pca_res$x) == input$comp1)
      comp2 <- which(colnames(pca_res$x) == input$comp2)
      
      # create PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca_plot <- ggplot(plot_tib, aes(x = PC1, y = PC2)) +
        geom_point() +
        labs(title = "Principal Component Analysis Plot",
             x = paste0(input$comp1, " (", round(var_explained[comp1], 2), "% variance)"),
             y = paste0( input$comp2, " (", round(var_explained[comp2], 2), "% variance)")) +
        theme_bw()
      return(pca_plot)
    }
  })
  
  observeEvent(input$submit_btn, {
    output$summary_data <- renderDataTable({
      if (!is.null(data_counts())) {
        processed_data <- summary_data(data_counts(), input$gvar, input$nz)
        processed_data
      }
    }, options = list(lengthChange = FALSE))
  })

}
#Render the whole thing
shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "Rshiny Project BF591"),
    sidebar,
    body
  ),
  server = server
)