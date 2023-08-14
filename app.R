library(shiny)
library(shinythemes)
library(DT)

rm(list=ls())

funcs_path <- "G:/Shared drives/Denby Lab Team Drive/Lab members/Harry/lettuce_data_shiny_app/"
source(file.path(funcs_path, 'sql_db/db_connect.R'))
source(file.path(funcs_path, 'plotting_funcs/plot_divset_cor.R'))
source(file.path(funcs_path, 'plotting_funcs/identify_subset_hubs.R'))

ui <- fluidPage(
  theme = shinytheme("cosmo"),  # Theme for the app

  titlePanel("Lettuce Data Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput("experiment", "Choose an Experiment",
                  choices = c("Lesion Size Correlation", "Time-Series Expression", "Gene Regulatory Network Analysis"),
                  selected = "Lesion Size Correlation"),

      selectInput("geneSelectionMethod", "Gene selection method",
                  choices = c("Lettuce GeneID", "Ortholog of Arabidopsis Genes", "Genes with GO-term", "Genes with Protein Domain"),
                  selected = "Genes with GO-term"),

      # Dynamic Default for Text Input based on Gene Selection Method
      uiOutput("geneInputUI"),

      # Dynamic UI elements based on experiment choice
      uiOutput("experimentFilters"),

      uiOutput("additionalFilters"),

      # Add actionButton to generate results
      actionButton("btn_generate", "Generate Results")
    ),


    mainPanel(
      # Display information if button is not clicked
      conditionalPanel(
        condition = "input.btn_generate == 0",
        tags$h3("Welcome to Lettuce Data Explorer!"),
        tags$p("This application allows you to perform integrated analyses on Lettuce data."),
        tags$p("Start by selecting an experiment and the gene selection method. You can then enter the specific genes or criteria you are interested in and adjust any additional filters or settings. Once you are ready, click the 'Generate Results' button to view the corresponding plots and data tables.")
      ),

      # Display the appropriate output based on the selected experiment after the button is clicked
      conditionalPanel(
        condition = "input.experiment == 'Lesion Size Correlation' &&  input.btn_generate != 0",
        plotOutput("mainPlot")
      ),
      # conditionalPanel(
      #   condition = "input.experiment == 'Lesion Size Correlation' && input.lesion_cor_plotType == 'Heatmap' && input.btn_generate != 0",
      #   plotOutput("heatmapPlot", height = "600px", width = "100%")
      # ),
      # conditionalPanel(
      #   condition = "input.experiment == 'Time-Series Expression' && input.btn_generate != 0",
      #   plotOutput("mainPlot")
      # ),
      conditionalPanel(
        condition = "input.experiment == 'Gene Regulatory Network Analysis' && input.btn_generate != 0",
        DTOutput("mainTable")
      )
    )
  ),

  # Footer with styled citations
  tags$footer(
    tags$hr(),
    tags$div(
      style = "border: 2px solid #ddd; padding: 10px; background-color: #f9f9f9; font-size: 16px;",  # Styles for the citation box
      tags$p(
        "Citations:",
        tags$br(),
        HTML("<a href='https://doi.org/10.1007/s00122-022-04129-5'>Pink, H., Talbot, A., Graceson, A. et al. Identification of genetic loci in lettuce mediating quantitative resistance to fungal pathogens. Theor Appl Genet, 135, 2481–2500 (2022).</a>"),
        tags$br(),
        HTML("<a href='https://doi.org/10.1101/2023.07.19.549542'>Pink, H., Talbot, A., Carter, R.,. et al. Identification of Lactuca sativa transcription factors impacting resistance to Botrytis cinerea through predictive network inference. bioRxiv (2023).</a>")
      )
    )
  )
)


server <- function(input, output, session) {

  # Render the appropriate default text based on gene selection method
  output$geneInputUI <- renderUI({
    switch(input$geneSelectionMethod,
           "Lettuce GeneID" = textInput("geneInput", "Enter lettuce GeneIDs as comma separated values",
                                        value="Lsat_1_v5_gn_7_33721, Lsat_1_v5_gn_8_115261, Lsat_1_v5_gn_8_116421, Lsat_1_v5_gn_8_116340"),
           "Ortholog of Arabidopsis Genes" = textInput("geneInput", "Enter Arabidopsis gene names or IDs as comma separated values",
                                                       value="ERF1, WRKY33, AT1G55020, AT5G42650"),
           "Genes with GO-term" = textInput("geneInput", "Enter GO Term",
                                            value="jasmonic acid mediated signaling pathway"),
           "Genes with Protein Domain" = textInput("geneInput", "Enter Protein Domain Description or Pfam ID", value="Cytochrome P450")
    )
  })

  input_vector <- reactive({
    if(!is.null(input$geneInput)) {
      # Split the input based on commas and whitespace
      genes <- unlist(strsplit(input$geneInput, split = ",\\s*"))
      return(genes)
    } else {
      return(NULL)
    }
  })

  input_list <- reactive({
    list(
      GeneIDs = if(input$geneSelectionMethod == 'Lettuce GeneID') input_vector() else NULL,
      At_orthologs = if(input$geneSelectionMethod == 'Ortholog of Arabidopsis Genes') input_vector() else NULL,
      GO_id = if(input$geneSelectionMethod == 'Genes with GO-term') input_vector() else NULL,
      protein_domain = if(input$geneSelectionMethod == 'Genes with Protein Domain') input_vector() else NULL
    )
  })


  output$experimentFilters <- renderUI({
    if (input$experiment == "Lesion Size Correlation") {

      tagList(
        selectInput('lesion_cor_plotType', "Plot style",
                    choices = c('Heatmap', 'Single Panel Scatterplot', 'Multi-panel Scatterplot'),
                    selected = 'Heatmap'),
        selectInput('lesion_cor_fungi', "Plot lesion size - expression correlation in response to which fungi?",
                    choices = c('S. sclerotiorum only', 'B. cinerea AND S. sclerotiorum'),
                    selected = 'S. sclerotiorum only'),
        sliderInput("topn_by_lesion_cor",
                    "Maximum of number of genes to return",
                    min = 1,
                    max = 100,
                    value = 10),
        sliderInput('lesion_corr_facet_title_size',
                    "Plot options: Gene label size",
                    min = 7,
                    max = 20,
                    value = 10)
      )
      # } else if (input$experiment == "Time-Series Expression") {
      #   selectInput("geneFilter2", "Filter by Time", choices = c("Time1", "Time2"))
    } else if (input$experiment == "Gene Regulatory Network Analysis") {
      tagList(
        selectInput('grn_output_type','Network table output',
                    choices = c('Individual edges','Aggregated regulator statistics'),
                    selected = 'Aggregated regulator statistics'),
        sliderInput("n_TFs",
                    "Maximum of Predicted Regulators to return:",
                    min = 1,
                    max = 30,
                    value = 10),
        sliderInput("min_targets",
                    "Minimum number of predicted targets within gene subset",
                    min = 1,
                    max = 50,
                    value = 8)
      )
    }
  })


  output$additionalFilters <- renderUI({
    if(input$lesion_cor_plotType =='Multi-panel Scatterplot'){
      tagList(
        sliderInput("lesion_corr_facet_nrows",
                    "Plot options: Number of grid rows",
                    min = 1,
                    max = 8,
                    value = 2),
        selectInput('lesion_corr_facet_scales',
                    'Plot options: Y-axis scale',
                    choices = c('Consistent Y-axis across all genes', 'Gene-specific Y-axis'),
                    selected = 'Consistent Y-axis across all genes')
      )

    }else if(input$grn_output_type == 'Aggregated regulator statistics'){
      tagList(
        selectInput("hubs_sortby", "Sort Predicted Regulators by;",
                    choices = c("FoldEnrichment","SubsetOutdegrees", "SumImportance"),
                    selected = "SumImportance")
      )

    }else{
      NULL
    }
  })


  ##initialize value for the plot

  make_plot <- reactiveVal()

  observeEvent(input$btn_generate, {
    # Update the table data when button is clicked
    # Check if the button was clicked
    if (input$experiment == "Lesion Size Correlation") {
      if(input$lesion_cor_fungi=='S. sclerotiorum only'){
        input_fungi = c('Scl')
      }else{
        input_fungi <- c('Scl','Bot')
      }
      args <- c(input_list(),
                list(fungi = input_fungi,
                     top_n_by_lesion_cor = as.integer(input$topn_by_lesion_cor),
                     return_heatmap = ifelse(input$lesion_cor_plotType =='Heatmap', TRUE, FALSE),
                     facet_title_size=input$lesion_corr_facet_title_size)
      )
      if(input$lesion_cor_plotType =='Multi-panel Scatterplot'){
        print('multi-panel scatter selected')
        args <- c(args,
                  list(single_panel=FALSE,
                       facet_rows=input$lesion_corr_facet_nrows,
                       facet_scales = ifelse(input$lesion_corr_facet_scales=='Gene-specific Y-axis',
                                             'free', 'free_x')))
      }

      p <- do.call(plot_divset_cor, args)
      }else if (input$experiment == "Time-Series Expression") {
            plot(sin(seq(1, 10, by = 0.1)), type = "l", main="Time-Series Expression")
      }else{
        p <- ggplot()
      }

    make_plot(p) ##store plot in the reactive value
    })




    output$mainPlot <- renderPlot({
      make_plot()
    }, height = 600, width = function() { return(session$clientData$output_mainPlot_width) })
  #
  #   # Initialize a reactiveValue for the table data
  tableData <- reactiveVal()
  # #
  observeEvent(input$btn_generate, {
    # Update the table data when button is clicked
    if (input$experiment == "Gene Regulatory Network Analysis") {
      if(trimws(input$geneInput)==""){
        df <- db_query("SELECT TF, n_targets as TotalOutdegrees,
                        CASE
                       WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN TF
                       ELSE 'Ls' || At_ShortName
                       END AS name
                       FROM grn_hubs WHERE TotalOutdegrees>=40")
      }else{
        args <- c(input_list(),
                  list(n_TFs = input$n_TFs,
                       min_subset_targets = input$min_targets,
                       order_by_column =  ifelse(input$grn_output_type=='Individual edges','SumImportance',input$hubs_sortby),
                       return_edges = ifelse(input$grn_output_type=='Individual edges',TRUE,FALSE)))
        df <- do.call(get_pred_regulators, args)
      }

      tableData(df) # Store the data frame in the reactiveVal
    }
  })

  # Render the table based on tableData reactiveVal
  output$mainTable <- renderDT({
    datatable(tableData())
  })

}


shinyApp(ui = ui, server = server)
