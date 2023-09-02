library(shiny)
library(shinythemes)
library(shinycssloaders)
library(DT)

rm(list=ls())


source('sql_db/db_connect.R')
source('plotting_funcs/helper_funcs.R')
source('plotting_funcs/plot_divset_cor.R')
source('plotting_funcs/identify_subset_hubs.R')
source('plotting_funcs/plot_timeseries_expr.R')

troubleshoot_prints <- FALSE

ui <- fluidPage(
  theme = shinytheme("cosmo"),  # Theme for the app

  titlePanel("Lettuce Data Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput("experiment", "Choose an Experiment",
                  choices = c( "Time-Series Expression", "Gene Regulatory Network Analysis","Lesion Size Correlation"),
                  selected = "Time-Series Expression"),

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

        tags$div(style = "font-size: 18px; line-height: 1.6;",   # This applies the styling to the entire div content

                 tags$h3("Welcome to Lettuce Data Explorer!"),
                 tags$p("This tool simplifies the exploration of our lettuce transcriptomic datasets. Quickly find lettuce genes using Arabidopsis symbols/IDs, GO terms, or protein domains, allowing you to focus on your genes of interest."),

                 tags$h4("Instructions:"),

                 tags$ul(
                   tags$li("Select an experiment to analyse Lettuce data:",

                           tags$ul(
                             tags$li(
                               strong("Time-series expression:"),
                               "Monitor lettuce gene responses to ", HTML("<i>Botrytis cinerea</i>"), " and ", HTML("<i>Sclerotinia sclerotioruim</i>"), ". View gene expression profiles."
                             ),
                             tags$li(
                               strong("Gene Regulatory Network Analysis:"),
                               "Explore transcriptional shifts during necrotrophic infections. View regulators or edges for selected genes."
                             ),
                             tags$li(
                               strong("Lesion Size Correlation:"),
                               "Examine gene expression's correlation with susceptibility (lesion size) to pathogens."
                             )
                           )),

                   tags$li("Choose a gene selection method:",
                           tags$ul(
                             tags$li("Lettuce GeneID: Input a Lettuce gene ID."),
                             tags$li("Ortholog of Arabidopsis Genes: Provide an Arabidopsis gene ID to find Lettuce ortholog."),
                             tags$li("Genes with GO-term: Enter a GO term for associated Lettuce genes."),
                             tags$li("Genes with Protein Domain: Search genes by protein domains.")
                           )),
                   tags$li("Specify genes or criteria based on your selected gene method."),
                   tags$li("Adjust filters or settings for your selected experiment."),
                   tags$li("Click 'Generate Results' to view plots and data tables.")
                 )
        )
      )

      ,

      # Display the appropriate output based on the selected experiment after the button is clicked
      # Lesion Size Correlation
      fluidRow(
        column(6,
               conditionalPanel(
                 condition = "(input.experiment == 'Lesion Size Correlation' || input.experiment == 'Time-Series Expression') && input.btn_generate != 0",
                 downloadButton('downloadData', 'Download Expression Data')
               )
        ),
        column(6,
               conditionalPanel(
                 condition = "input.experiment == 'Lesion Size Correlation' &&  input.btn_generate != 0",
                 downloadButton('downloadLesionSizeData', 'Download Lesion Size Data')
               )
        )
      ),
      conditionalPanel(
        condition = "input.experiment == 'Lesion Size Correlation' &&  input.btn_generate != 0",
        withSpinner(plotOutput("lesionCorrPlot"))
      ),
      # conditionalPanel(
      #   condition = "input.experiment == 'Lesion Size Correlation' && input.btn_generate != 0",
      #   withSpinner(DTOutput("lesionCorrTable"))
      # ),

      # Time-Series Expression
      conditionalPanel(
        condition = "input.experiment == 'Time-Series Expression' && input.btn_generate != 0",
        withSpinner(plotOutput("timeSeriesPlot"))
      ),
      # conditionalPanel(
      #   condition = "input.experiment == 'Time-Series Expression' && input.btn_generate != 0",
      #   withSpinner(DTOutput("timeSeriesTable"))
      # ),

      # Gene Regulatory Network Analysis
      conditionalPanel(
        condition = "input.experiment == 'Gene Regulatory Network Analysis' && input.btn_generate != 0",
        withSpinner(DTOutput("mainTable"))
      )
    )
  ),  # This closes the sidebarLayout

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
           "Genes with GO-term" = textInput("geneInput", "Enter GO Term Description or ID",
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
        selectInput("signif_cor_only", "Do you wish to return significant lesion size correlated genes only?",
                    choices = c('Yes'=TRUE,'No'=FALSE),
                    selected = TRUE),
        selectInput("DEGs_filter",
                    "Do you wish to remove genes with confounding differential expression? See Pink et al., (2022) for more details",
                    choices = c("Yes (Recommended)" = TRUE, "No" = FALSE),
                    selected = TRUE),
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
      } else if (input$experiment == "Time-Series Expression") {
        tagList(
          selectInput('degs_input',
                      'Do you wish to filter gene based on differenetial expression',
                      choices = c("No, show all genes","Only display DEGs",
                                  "Only display upregulated genes","Only display downregulated genes"),
                      selected = "Only display DEGs"),
          selectInput('timeseries_fungi', "Plot time-series expression in response to which fungi?",
                      choices = c('S. sclerotiorum only','B. cinerea only' ,'B. cinerea AND S. sclerotiorum'),
                      selected = 'B. cinerea AND S. sclerotiorum'),
          numericInput("max_n",
                      "Maximum of number of genes to plot (filtered by MaxAbsLog2FC)",
                      min = 1,
                      max = 1000,
                      value = 1000),
          selectInput('time_series_plot_type', 'Plot options: Plot style',
                      choices = c('Heatmap', 'Single Panel Line Plot', 'Multi Panel Line Plot'),
                      selected = 'Heatmap'),
          selectInput('include_mock','Plot options: Do you wish to include mock gene expression',
                      choices = c('Yes','No'),
                      selected = 'Yes'),
          sliderInput('time_series_label_size','Plot option: gene label size',
                      min = 7,
                      max = 20,
                      value = 10)

        )
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
  if(input$experiment == "Lesion Size Correlation"){
    if(isTruthy(input$lesion_cor_plotType) ){
      if(input$lesion_cor_plotType =='Multi-panel Scatterplot'){
        tagList(
          sliderInput("lesion_corr_facet_nrows",
                      "Plot options: Number of grid rows",
                      min = 1,
                      max = 8,
                      value = 2),
          selectInput('lesion_corr_facet_scales',
                      'Plot options: Y-axis scale',
                      choices = c('Consistent Y-axis across all genes'='free_x',
                                  'Gene-specific Y-axis'='free'),
                      selected = 'free_x')
        )}
    }}else if(input$experiment=="Time-Series Expression"){
      if(isTruthy(input$time_series_plot_type)){
        if(input$time_series_plot_type =='Multi Panel Line Plot'){
          tagList(
            sliderInput("time_series_facet_nrows",
                        "Plot options: Number of grid rows",
                        min = 1,
                        max = 8,
                        value = 2),
            selectInput('timeseries_facet_scales',
                        'Plot options: Y-axis scale',
                        choices = c('Consistent Y-axis across all genes', 'Gene-specific Y-axis'),
                        selected = 'Consistent Y-axis across all genes')
          )}
        }
    }else if(input$experiment == "Gene Regulatory Network Analysis"){


      if(troubleshoot_prints){ print('recognised valid grn_output_type')}

      if(isTruthy(input$grn_output_type)){

        if(input$grn_output_type == 'Aggregated regulator statistics'){
          tagList(
            selectInput("hubs_sortby", "Sort Predicted Regulators by;",
                        choices = c("FoldEnrichment","SubsetOutdegrees", "SumImportance"),
                        selected = "SumImportance")
          )
        }

      }



    }else{
      if(troubleshoot_prints){ print('no ')}
      NULL
    }
  })


  ##initialize value for the plot

  make_plot <- reactiveVal()
  plot_data <- reactiveVal()
  lesion_size_data <- reactiveVal()

  observeEvent(input$btn_generate, {
    # Update the table data when button is clicked
    # Check if the button was clicked

    if(troubleshoot_prints){
      print(input$experiment)
      print(input$geneInput)
      print(input$geneSelectionMethod)
      print(input$lesion_cor_plotType)
      print(input$lesion_cor_fungi)
      print(input$grn_output_type)
    }
    p<- ggplot()
    df <- data.frame()
    divset_lesion_data <- data.frame()



    ##initiate args with input list

    if (input$experiment == "Lesion Size Correlation") {
      if(isTruthy(input$lesion_cor_fungi)){
        if(troubleshoot_prints){ print('recognised lesion cor fungi')}
        if(input$lesion_cor_fungi=='S. sclerotiorum only'){
          input_fungi = c('Scl')
          }
        }else{
          input_fungi <- c('Scl','Bot')

        }


      if(troubleshoot_prints){ print(input_fungi)}



      args <- c(input_list(),
                list(fungi = input_fungi,
                     top_n_by_lesion_cor = as.integer(input$topn_by_lesion_cor),
                     signif_cor_only = input$signif_cor_only,
                     DEGs_filtered= input$DEGs_filter,
                     return_heatmap = ifelse(input$lesion_cor_plotType =='Heatmap', TRUE, FALSE),
                     facet_title_size=input$lesion_corr_facet_title_size) )
      if(isTruthy(input$lesion_cor_plotType)){
        if(input$lesion_cor_plotType =='Multi-panel Scatterplot'){
        args <- c(args,
                  list(single_panel=FALSE,
                       facet_rows=input$lesion_corr_facet_nrows,
                       facet_scales = input$lesion_corr_facet_scales))
        }
        }

      if(length(args)>0){
        divset_plot_obj <- do.call(plot_divset_cor, args)
        p <- divset_plot_obj$plot
        df <- divset_plot_obj$exp_data
        divset_lesion_data <- divset_plot_obj$lesion_data
        if(troubleshoot_prints){ print(args)}
      }


      }else if(input$experiment == "Time-Series Expression"){

        ##get fungi input vector
        if(isTruthy(input$timeseries_fungi)){
          input_fungi <- c()
          if(grepl('S. scl',input$timeseries_fungi)){
            input_fungi <- c(input_fungi, 'sclero')
          }
          if(grepl('B. cin',input$timeseries_fungi)){
            input_fungi <- c(input_fungi, 'bot')
          }
        }else{
          input_fungi <- c('bot','sclero')
        }

        ##get degs input
        up_down_input = 'both'
        if(isTruthy(input$degs_input)){
          if(grepl('upregulated',input$degs_input)){
            up_down_input = 'up'
          }else if(grepl('downregulated',input$degs_input)){
            up_down_input = 'down'
          }
        }

        input_plot_type <- 'line'
        single_panel = FALSE
        input_facet_nrow =2
        input_facet_scales = 'free_y'
        if(isTruthy(input$time_series_plot_type)){
          if(input$time_series_plot_type =="Heatmap"){
            input_plot_type = 'heatmap'
          }else if(input$time_series_plot_type=="Single Panel Line Plot"){
            single_panel = TRUE
          }else if(input$time_series_plot_type=="Multi Panel Line Plot"){
            input_facet_nrow = as.integer(input$time_series_facet_nrows)
            input_facet_scales = ifelse(input$timeseries_facet_scales=='Gene-specific Y-axis',
                                  'free_y', 'free_x')
          }

        }

        if(troubleshoot_prints){ print('setting time-series args')}



        args <- c(input_list(),
                  list(max_n = as.integer(input$max_n),
                       fungi = input_fungi,
                       include_mock = input$include_mock =='Yes',
                       overlap_DEGs_only = !input$degs_input =="No, show all genes",
                       up_down_only = up_down_input,
                       plot_type = input_plot_type,
                       strip.text.x_size = input$time_series_label_size,
                       single_panel = single_panel,
                       facet_nrow=input_facet_nrow,
                       facet_scales=input_facet_scales))


        if(length(args)>0){
          timeseries_plot_obj <- do.call(plot_timeseries_expr, args)
          df <- timeseries_plot_obj$data
          p <- timeseries_plot_obj$plot

          if(troubleshoot_prints){ print(args)}
        }


      }

    if(troubleshoot_prints){ print('assigning plot to reactive variable')}
    make_plot(p) ##store plot in the reactive value
    plot_data(df)
    lesion_size_data(divset_lesion_data)
    if(troubleshoot_prints){ print('plot assigned')}
    })




  output$lesionCorrPlot <- renderPlot({
    make_plot()
  }, height = 710, width = function() { return(session$clientData$output_lesionCorrPlot_width) })



  output$timeSeriesPlot <- renderPlot({
    make_plot()
  },
  height = 710,
  width = function() { return(session$clientData$output_timeSeriesPlot_width) })


  output$downloadData <- downloadHandler(
    filename = function() {
      paste(gsub(' ','-',tolower(input$experiment)),'-expression-data-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(plot_data(), file)
    }
  )
  # New download handler specifically for Lesion Size Data
  output$downloadLesionSizeData <- downloadHandler(
    filename = function() {
      paste('lesion-size-data-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(lesion_size_data(), file)
    }
  )



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
