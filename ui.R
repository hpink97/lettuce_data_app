# Import necessary libraries
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)

ui <- fluidPage(

  # Set the theme for the app
  theme = shinytheme("cosmo"),

  # Title of the web page
  titlePanel("Lettuce Data Explorer"),

  # Add a div section for citations banner at the top of the app
  tags$div(
    # Styling for the citation box
    style = "border: 2px solid #ddd; padding: 10px; background-color: #f9f9f9; font-size: 16px;",
    tags$p(
      "Citations:",  # Title for the citations
      tags$br(),  # Break line
      # Link to the first citation
      HTML("<a href='https://doi.org/10.1007/s00122-022-04129-5'>Pink, H., Talbot, A., Graceson, A. et al. Identification of genetic loci in lettuce mediating quantitative resistance to fungal pathogens. Theor Appl Genet, 135, 2481–2500 (2022).</a>"),
      tags$br(),  # Break line
      # Link to the second citation
      HTML("<a href='https://doi.org/10.1101/2023.07.19.549542'>Pink, H., Talbot, A., Carter, R.,. et al. Identification of Lactuca sativa transcription factors impacting resistance to Botrytis cinerea through predictive network inference. bioRxiv (2023).</a>")
    )
  ),

  # Structure layout to have a sidebar and a main content panel
  sidebarLayout(

    # Define the sidebar contents
    sidebarPanel(
      # Dropdown for selecting dataset for analysis
      selectInput("experiment", "Select a dataset to analyse",
                  choices = c( "Time-Series Expression", "Gene Regulatory Network Analysis","Lesion Size Correlation")),

      # Dropdown for gene selection method
      selectInput("geneSelectionMethod", "Gene selection method",
                  choices = c("Lettuce GeneID", "Ortholog of Arabidopsis Genes", "Genes with GO-term", "Genes with Protein Domain")),

      # Further dynamic user inputs that changes based on previous selection
      # Dropdowns specific to gene selection method
      uiOutput("geneInputUI"),

      # Further dynamic dropdowns based on experiment/dataset selected
      uiOutput("experimentFilters"),

      # Further dynamic dropdowns based on experimentFilters selected
      uiOutput("additionalFilters"),

      # Button to trigger generation of results
      actionButton("btn_generate", "Generate Results")
    ),

    # Define the main content panel
    mainPanel(

      # Display welcome information if button is not clicked
      conditionalPanel(
        condition = "input.btn_generate == 0",
        tags$div(style = "font-size: 18px; line-height: 1.6;",
                 # Welcome and brief introduction
                 tags$h3("Welcome to Lettuce Data Explorer!", style="font-weight: bold;"),
                 tags$p("This tool simplifies the exploration of our lettuce transcriptomic datasets. Quickly find lettuce genes using Arabidopsis symbols/IDs, GO terms, or protein domains, allowing you to focus on your genes of interest."),

                 # Detailed steps and descriptions
                 tags$h4("Step 1: Select a dataset to analyse", style="font-weight: bold;"),
                 tags$ul(
                   tags$li(strong("Time-series expression:"),
                           "Explore dynamic expression profiles ..."),
                   tags$li(strong("Gene Regulatory Network Analysis:"),
                           "Explore predicted transcriptional regulation ..."),
                   tags$li(strong("Lesion Size Correlation:"), "Examine gene expression's correlation ...")
                 ),

                 tags$h4("Step 2: Select your genes of interest ...", style="font-weight: bold;"),
                 tags$ul(
                   tags$li(strong("Lettuce GeneID:"), "Input Lettuce gene ID(s)."),
                   tags$li(strong("Orthologues of Arabidopsis Genes:"), "Provide an Arabidopsis gene ID ..."),
                   tags$li(strong("Genes with GO-term:"), "Enter a GO term ..."),
                   tags$li(strong("Genes with Protein Domain:"), "Search genes ...")
                 ),
                 tags$h4("Step 3: Dataset-specific options", style="font-weight: bold;"),
                 tags$ul(
                   tags$li("Dataset-specific gene selection criteria ..."),
                   tags$li("Plot customisation options")
                 ),

                 tags$h4("Step 4: Click 'Generate Results'", style="font-weight: bold;")
        )
      ),

      # Data and plot download buttons
      fluidRow(
        column(4,
               ## Show download expression button if experiment chosen is Lesion correlation or time-series expression
               conditionalPanel(
                 condition = "(input.experiment == 'Lesion Size Correlation' || input.experiment == 'Time-Series Expression') && input.btn_generate != 0",
                 downloadButton('downloadData', 'Download Expression Data')
               )
        ),
        column(4,
               ## Show download plot button if experiment chosen is Lesion correlation or time-series expression
               conditionalPanel(
                 condition = "(input.experiment == 'Lesion Size Correlation' || input.experiment == 'Time-Series Expression') && input.btn_generate != 0",
                 actionButton("showDownloadOptions", "Download Plot")
               )
        ),
        column(4,
               ## Show download lesion correlation button if experiment chosen is Lesion correlation
               conditionalPanel(
                 condition = "input.experiment == 'Lesion Size Correlation' &&  input.btn_generate != 0",
                 downloadButton('downloadLesionSizeData', 'Download Lesion Size Data')
               )
        )
      ),

      # Conditional displays of plots and tables
      conditionalPanel(
        condition = "input.experiment == 'Lesion Size Correlation' &&  input.btn_generate != 0",
        withSpinner(plotOutput("lesionCorrPlot"))
      ),
      conditionalPanel(
        condition = "input.experiment == 'Time-Series Expression' && input.btn_generate != 0",
        withSpinner(plotOutput("timeSeriesPlot"))
      ),
      conditionalPanel(
        condition = "input.experiment == 'Gene Regulatory Network Analysis' && input.btn_generate != 0",
        withSpinner(DTOutput("mainTable")),
        uiOutput("dynamicLegend")
      )
    )
  )  # This closes the sidebarLayout
)
