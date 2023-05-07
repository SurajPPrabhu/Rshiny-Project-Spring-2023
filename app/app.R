library(shiny)
library(shinydashboard)
library(tidyverse)
library(scales)
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
                                     min = -50, max = 100, step = 1, value = 50),
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
      read_table(input$file2$datapath)
    }
  })
  
  # Define summary_data function to process counts data
  summary_data <- function(data_counts,gvar, nz) {
    # Filter genes based on variance and number of non-zero samples
    genes_var_filtered <- data_counts[rowSums(data_counts > 0) >= nz & apply(data_counts, 1, var) >= gvar, ]
    
    # Get summary statistics for the filtered data
    num_rows_filtered <- nrow(genes_var_filtered)
    num_rows_not_filtered <- nrow(data_counts)-num_rows_filtered
    pgene <- percent((num_rows_filtered/nrow(data_counts)))
    fgene<-percent((num_rows_not_filtered)/nrow(data_counts))
    
    # Create summary data frame
    summary_df <- data.frame(Description = c("Number of Genes", "Number of Samples", "Number of Genes(Filtered):","% of genes filtered", "Number of Genes not filtered:","%genes not filtered:"), Value = c(nrow(data_counts), ncol(data_counts), num_rows_filtered,pgene, num_rows_not_filtered,fgene))
    
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
        tabPanel("Counts Information"),
        tabPanel("Counts Plot")
      )
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