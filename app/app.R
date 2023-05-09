library(shiny)
library(shinydashboard)
library(tidyverse)
library(scales)
library(gplots)
library(RColorBrewer)
library(shinyWidgets)
library(colourpicker)
library(gridExtra)
library(glue)
library(shinythemes)
library(DT)
library(ggbeeswarm)

#defining the layout and sidebars
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Sample Information Exploration", tabName = "data"),
    menuItem("Counts", tabName = "counts"),
    menuItem("Differential Expression", tabName ="diffexp" ),
    menuItem("Individual gene Analysis", tabName = "ind_gene")
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
            fluidPage(sidebarPanel(
              fileInput("file2", "Upload CSV file containing counts information"),
              sliderInput("gvar", "Select the minimum percentile of variance desired in the genes",
                          min = -50, max = 100, value = 50),
              sliderInput("nz", "Select the minimum number of nonzero samples",
                          min = 0, max = 100, value = 50),
              actionButton("submit_btn", "Submit")
            ),
            mainPanel(
              uiOutput("counts_tabs")
            ))
    )
    ,
    tabItem(tabName = "diffexp", 
            sidebarLayout(
              sidebarPanel(
                fileInput(inputId = "file3","Load differential expression results:", placeholder = "deseq_res.csv", accept = ".csv"), 
                radioButtons(inputId = "x_var", "Choose the column for x-axis", choices=c("baseMean","log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                radioButtons(inputId = "y_var", "Choose the column for y-axis", choices=c("baseMean","log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
                colourInput(inputId = "color1", "Select base point color", value = "midnightblue", showColour = c("both", "text", "background"),
                            palette = c("square", "limited"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE),
                colourInput(inputId = "color2", "Highlight point color", value = "springgreen", showColour = c("both", "text", "background"), 
                            palette = c("square", "limited"), allowedCols = NULL, allowTransparent = FALSE, returnName = FALSE, closeOnClick = FALSE),
                sliderInput(inputId = "slider_value", "Select the magnitude of the p adjusted coloring:", min = -100, max = 0, value = -50),
                actionButton("plot", "Submit")
              ),
              mainPanel(
                h3("Differential Expression Analysis"),
                p("Load differential expression analysis results for visualization.")
              )
            ),uiOutput("diff_exp")
    ),
    
    tabItem(tabName = "ind_gene",
            sidebarLayout(
              sidebarPanel(fileInput(inputId = 'counts_id', label = 'Load normalized counts matrix CSV'),
                           fileInput(inputId = 'meta_id', label = 'Load sample information matrix CSV'),
                           uiOutput("meta_selector"),
                           #insert gene search box here
                           textInput("gene", label = "Enter gene to search for", placeholder = "ENSG00000000003.10"),
                           selectInput("plotType", label = "Choose what type of plot to make", choices = c("Bar", "Violin", "Beeswarm")),
                           actionBttn("go", "Submit")
              ),
              mainPanel(
                plotOutput("various_plot")
                
              )
            )
            
    ))
  
)



#defining server side functions
server <- function(input, output) {
  options(shiny.maxRequestSize = 30*1024^2)
  
  # Read metadata file
  data_summary <- reactive({
    if (is.null(input$file1)) {
      return(NULL)
    } else {
      read_csv(input$file1$datapath)
    }
  })
  
  # Define data_counts reactive function to read counts file
  data_counts <- reactive({
    if (is.null(input$file2)) {
      return(NULL)
    } else {
      read_csv(input$file2$datapath, col_names = TRUE)
    }
  })
  
  # Function to load diff data 
  diff_data <- reactive({
    if (is.null(input$file3)) {
      return(NULL)
    } else {
      read_csv(input$file3$datapath)
    }
  })
  
  data_sum <- reactive({
    if (is.null(input$file4)) {
      return(NULL)
    } else {
      read_csv(input$file4$datapath)
    }
  })
  
  data_ct <- reactive({
    if (is.null(input$file5)) {
      return(NULL)
    } else {
      read_csv(input$file5$datapath)
    }
  })
  
  #Define summary tab for the meta data
  meta_summary <- function(data) {
    col_types <- sapply(data, class)
    col_means <- sapply(data, mean, na.rm = TRUE)
    col_sds <- sapply(data, sd, na.rm = TRUE)
    summ_tib <- tibble("Columns" = names(data), "Type" = col_types, 
                       "Mean" = col_means, "SD(+/-)" = col_sds)
    return(summ_tib)
  }
  
  
  
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
  # summary_data <- function(data_counts, perc_var, nz_genes) {
  #   # Calculate variance percentile threshold
  #   var_threshold <- quantile(apply(data_counts, 1, var, na.rm='TRUE'), probs = perc_var/100)
  #   
  #   # Filter genes based on variance and number of non-zero samples
  #   genes_var_filtered <- data_counts[rowSums(data_counts > 0) >= nz_genes & apply(data_counts, 1, var) >= var_threshold, ]
  #   
  #   # Get summary statistics for the filtered data
  #   num_rows_filtered <- nrow(genes_var_filtered)
  #   num_rows_not_filtered <- nrow(data_counts) - num_rows_filtered
  #   pgene <- percent((num_rows_filtered / nrow(data_counts)))
  #   fgene <- percent((num_rows_not_filtered) / nrow(data_counts))
  #   
  #   # Create summary data frame
  #   summary_df <- data.frame(Description = c("Number of Genes", "Number of Samples", "Number of Genes (Filtered):", "% of genes filtered", "Number of Genes not filtered:", "%genes not filtered:"), Value = c(nrow(data_counts), ncol(data_counts)-1, num_rows_filtered, pgene, num_rows_not_filtered, fgene))
  #   
  #   return(summary_df)
  # } 
  
  
  # Render the data tabs UI
  output$data_tabs <- renderUI({
    if (is.null(input$file1)) {
      return(NULL)
    } else {
      tabsetPanel(
        tabPanel("Sample Summary", tableOutput("summary")),
        tabPanel("View my Data", dataTableOutput("viewall")),
        tabPanel("Plot my data", plotOutput("density_plot"))
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
        tabPanel("PCA" ,selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),selected="PC1"),
                 selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2")
                 , plotOutput("plot_pca"))
      )
    }
  })
  
  output$diff_exp<-renderUI({
    if (is.null(diff_data())) {
      return(NULL)
    } else {tabsetPanel(tabPanel("Table",dataTableOutput("deg_table")),
                        tabPanel("Volano Plot", plotOutput("volcano")),
                        tabPanel("Plot Table", dataTableOutput("v_table"))
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
  
  # Function to calculate PCA and create PCA plot
  calc_pca_plot <- function(data, gvar, comp1, comp2) {
    filt_tib <- data %>% 
      mutate(variance = apply(data[-1], MARGIN = 1, FUN = var), .after = gene) # calculate variance for filtering
    perc_val <- quantile(filt_tib$variance, probs = gvar/100, na.rm = TRUE)   # calculate percentile
    filt_tib <- filter(filt_tib, variance >= perc_val) # filter the tibble
    pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) # transpose the data and perform PCA
    
    # extract variance explained by each component
    var_explained <- summary(pca_res)$importance[2,] * 100
    
    # get selected principal components
    comp1 <- which(colnames(pca_res$x) == comp1)
    comp2 <- which(colnames(pca_res$x) == comp2)
    # get selected principal components
    
    
    # create PCA plot
    plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
    pca_plot <- ggplot(plot_tib, aes(x = PC1, y = PC2)) +
      geom_point() +
      labs(title = "Principal Component Analysis Plot",
           x = paste0(comp1, " (", round(var_explained[comp1], 2), "% variance)"),
           y = paste0(comp2, " (", round(var_explained[comp2], 2), "% variance)")) +
      theme_bw()
    return(pca_plot)
  }
  
  # Output render for PCA plot
  output$plot_pca <- renderPlot({
    if (is.null(data_counts())) {
      return(NULL)
    } else {
      pca_plot <- calc_pca_plot(data_counts(), input$gvar, input$comp1, input$comp2)
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
  
  density_plots <- function(data) {
    # Get index of numeric columns
    num_cols <- sapply(data, is.numeric)
    num_cols <- which(num_cols)
    
    if (length(num_cols) == 0) {
      return(NULL)
    }
    
    # Subset data frame to numeric columns
    data_num <- data[, num_cols]
    
    # Create density plots
    plots <- list()
    for (i in 1:length(num_cols)) {
      plots[[i]] <- ggplot(data = data_num, aes_string(x = paste0("`", names(data_num)[i], "`"))) + 
        geom_density() +
        theme_grey() +
        # labs(title = names(data_num)[i])
        labs(title = paste0("`", names(data_num)[i], "`"))
    }
    return(plots)}
  
  
  
  # Function to plot volcano plot using input filters from sliders 
  volc_plot<- function(data, x_var, y_var, slider_value, color1, color2) {
    data_filtered <- data %>% drop_na()
    volcano <- ggplot(data_filtered, aes(x= !!sym(x_var), y= -log10(!!sym(y_var)), color= color)) + 
      geom_point(aes(color = if_else(data_filtered[[y_var]] < 10^slider_value, TRUE, FALSE))) +   
      scale_color_manual(name = glue("{x_var} < 1*10^{slider_value}"),
                         values = c("FALSE" = color1,
                                    "TRUE" = color2),
                         labels = c("FALSE", "TRUE")) + 
      labs(title = "Volcano Plot", x = x_var, y = y_var) +
      theme_linedraw() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
    return(volcano)
  }
  
  # Function to generate table with data for DEGs filtered in volcano plot 
  volc_table <- function(file, slider) {
    # Filter data frame by adjusted p-value threshold
    new_df <- file[file$padj < 10^slider, ]
    
    # Format p-value columns to display more digits
    new_df$pvalue <- format(new_df$pvalue, digits = 5, scientific = TRUE)
    new_df$padj <- format(new_df$padj, digits = 5, scientific = TRUE)
    
    # Return filtered and formatted data frame
    return(new_df)
  }
  output$deg_table <- renderDataTable({
    diff_data()
  }, options = list(scrollX = TRUE))
  
  # Add observeEvent for submit button
  observeEvent(input$plot1, {
    # Generate plot and table with filtered data
    data3 <- diff_data()
    output$volcano <- renderPlot({
      volc_plot(data3, input$x_var, input$y_var, input$slider_value, input$color1, input$color2)
    })
    output$v_table <- renderDataTable({
      volc_table(data3, input$slider_value)
    }, options = list(scrollX = TRUE))
  })
  
  
  output$density_plot <- renderPlot({
    plots <- density_plots(data_summary())
    if (!is.null(plots)) {
      grid.arrange(grobs = plots, ncol = 2)
    }
  })
  
  output$summary<-renderTable(
    {meta_summary(data_summary())
      
    }
  )
  
  output$viewall<- renderDataTable({
    data_summary()
  }, options = list(scrollX = TRUE))
  
  
  ######4th tab
  
  # function to load normalized counts input file
  load_counts <- reactive({
    req(input$counts_id)
    read_csv(input$counts_id$datapath)
  })
  
  # function to load sample information input file
  load_meta <- reactive({
    req(input$meta_id)
    read_csv(input$meta_id$datapath)
  })
  
  # function to dynamically generate categorical variable options
  output$meta_selector <- renderUI({
    meta <- load_meta()
    if (!is.null(meta)) {
      cat_vars <- names(meta)[sapply(meta, is.numeric)]
      selectInput("metachoice", "Select categorical variable:", choices = cat_vars)
    }
  })
  # function to make distribution plots
  various_plots <- function(counts_tib, meta_tib, meta_cat, select_gene, plot_type) {
    counts_tib <- column_to_rownames(counts_tib, var = "gene")
    gene_counts <- as.numeric(as.vector(counts_tib[select_gene, ]))
    plot_tib <- tibble(Gene_Counts = gene_counts, meta_value = meta_tib[[meta_cat]])
    
    if (plot_type == "Beeswarm") {
      plot <- ggplot(plot_tib, aes(x = meta_value, y = Gene_Counts)) +
        geom_beeswarm() +
        theme_bw() +
        labs(title = str_c("Plot of gene counts vs ", meta_cat))
    } else if (plot_type == "Violin") {
      plot <- ggplot(plot_tib, aes(x = meta_value, y = Gene_Counts)) +
        geom_violin() +
        theme_bw() +
        labs(title = str_c("Plot of gene counts vs ", meta_cat))
    } else if (plot_type == "Bar") {
      plot <- ggplot(plot_tib, aes(x = meta_value)) +
        geom_bar() +
        theme_bw() +
        labs(title = str_c("Plot of gene counts vs ", meta_cat))
    }
    
    return(plot)
  }
  
  # methods to return outputs
  observeEvent(input$go, {
    output$various_plot <- renderPlot({
      various_plots(load_counts(), load_meta(), input$metachoice, input$gene, input$plotType)
    })
  })
}
#Render the whole thing
shinyApp(
  ui = fluidPage(dashboardPage(
    dashboardHeader(title = "Rshiny Project BF591"),
    sidebar,
    body
  )),
  server = server
)
