#load packages

library(shiny)
library(tidyverse)
library(DT)

# Set to TRUE if you need to print debug/troubleshoot messages
troubleshoot_prints <- FALSE


# Define server for shiny app
server <- function(input, output, session) {

  # Reactive value to store properties related to plot download
  plotDownloadUserChoices <- reactiveValues(fileType = NULL, imgWidth = NULL, imgHeight = NULL, imgResolution = NULL)

  # Initialize reactive values to store plot, plot data, lesion size data, and legend text

  make_plot <- reactiveVal(NULL)
  plot_data <- reactiveVal(NULL)
  lesion_size_data <- reactiveVal(NULL)
  tableData <- reactiveVal(NULL)

  LegendText <- reactiveVal(NULL)




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

  # Get the genes from input by splitting them on comma
  input_vector <- reactive({
    if(!is.null(input$geneInput)) {
      # Split the input based on commas and whitespace
      genes <- unlist(strsplit(input$geneInput, split = ",\\s*"))
      return(genes)
    } else {
      return(NULL)
    }
  })


  # Create a list of gene selection inputs
  # Selection methods not used will be null
  # Resulting list can be used directly in our plotting funcs

  input_list <- reactive({
    list(
      GeneIDs = if(input$geneSelectionMethod == 'Lettuce GeneID') input_vector() else NULL,
      At_orthologs = if(input$geneSelectionMethod == 'Ortholog of Arabidopsis Genes') input_vector() else NULL,
      GO_id = if(input$geneSelectionMethod == 'Genes with GO-term') input_vector() else NULL,
      protein_domain = if(input$geneSelectionMethod == 'Genes with Protein Domain') input_vector() else NULL
    )
  })


  ##########################################################################
  ########## Render UI that filters based on selected experiment  ##########
  ##########################################################################

  output$experimentFilters <- renderUI({

  ## Lesion size correlation specific filters
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


    ## Render Time-series expression specific UI filters
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

        sliderInput('time_series_label_size','Plot option: gene label size',
                    min = 7,
                    max = 20,
                    value = 10)

      )





    ## Render GRN-specific UI filters
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






  #################################################################################
  ###### Render further UI dropdowns/inputs based experiment specific inputs ######
  ###### Mainly plot options and table sorts                                 ######
  #################################################################################



  output$additionalFilters <- renderUI({
    if(input$experiment == "Lesion Size Correlation"){
      if(isTruthy(input$lesion_cor_plotType) ){

        ## Render UI specific to multi-panel scatterplots of Lesion-correlation
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

          #Render UI specific to multi-panel line plots of time-series expression
          if(input$time_series_plot_type =='Multi Panel Line Plot'){
            tagList(
              selectInput('include_mock','Plot options: Do you wish to include mock gene expression',
                          choices = c('Yes','No'),
                          selected = 'Yes'),
              sliderInput("time_series_facet_nrows",
                          "Plot options: Number of grid rows",
                          min = 1,
                          max = 8,
                          value = 2),
              selectInput('timeseries_facet_scales',
                          'Plot options: Y-axis scale',
                          choices = c('Consistent Y-axis across all genes', 'Gene-specific Y-axis'),
                          selected = 'Consistent Y-axis across all genes')
            )}else if(input$time_series_plot_type !='Single Panel Line Plot'){

          #Render UI specific to time-series expression plots that are neither single panel or ...
          #multi-panel line plots (i.e they are heatmaps)
              tagList(
                selectInput('include_mock','Plot options: Do you wish to include mock gene expression',
                            choices = c('Yes','No'),
                            selected = 'Yes'))
            }
        }
      }else if(input$experiment == "Gene Regulatory Network Analysis"){


        if(troubleshoot_prints){ print('recognised valid grn_output_type')}

        if(isTruthy(input$grn_output_type)){



          ## Render UI specific to GRN Hub gene data tables

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



  ##########################################################################
  ####  ACTIONS THAT TAKE PLACE WHEN GENERATE RESULTS BUTTON IS CLICKED ####
  ##########################################################################


  # The "observeEvent" function listens for an event, in this case, the event is a button click (input$btn_generate).
  observeEvent(input$btn_generate, {

    # If the "troubleshoot_prints" flag is TRUE, debugging input values.
    if(troubleshoot_prints){
      print(input$experiment)
      print(input$geneInput)
      print(input$geneSelectionMethod)
      print(input$lesion_cor_plotType)
      print(input$lesion_cor_fungi)
      print(input$grn_output_type)
    }

    ## initiate empty variables
    ## avoid errors filling reactive values with variables that don't exist
    p<- ggplot()
    df <- data.frame()
    grn_df <- data.frame
    divset_lesion_data <- data.frame()
    text <- ""






    # actions to be taken if experiment == lesion_cor when button clicked
    if (input$experiment == "Lesion Size Correlation") {


      ## Set fungal infection thta we require expression data for
      if(isTruthy(input$lesion_cor_fungi)){
        if(troubleshoot_prints){ print('recognised lesion cor fungi')}

        if(input$lesion_cor_fungi=='S. sclerotiorum only'){
          input_fungi = c('Scl')
        }
      }else{
        input_fungi <- c('Scl','Bot')

      }


      if(troubleshoot_prints){ print(input_fungi)}



      ##Create an argument list to provide to `plot_divset_cor` function
      ## fill with arguments that are consisnet across all plottypes

      args <- c(input_list(), #start with input containing gene selection criteria


                list(
                  ##user-selected input fungi vector
                  fungi = input_fungi,
                  ## how many genes to return
                  top_n_by_lesion_cor = as.integer(input$topn_by_lesion_cor),
                  ## boolean for whether to filter significant lesion cor only
                  signif_cor_only = input$signif_cor_only,
                  ## boolean for whether to filter for confounding DE
                  DEGs_filtered= input$DEGs_filter,
                  ## plot type - heatmap Y or N
                  return_heatmap = ifelse(input$lesion_cor_plotType =='Heatmap', TRUE, FALSE),
                  ## Plot output title size
                  facet_title_size=input$lesion_corr_facet_title_size)
                )

      ## Additional arguments specific to multi-panel scatterplot
      if(isTruthy(input$lesion_cor_plotType)){
        if(input$lesion_cor_plotType =='Multi-panel Scatterplot'){
          args <- c(args,
                    list(single_panel=FALSE,
                         facet_rows=input$lesion_corr_facet_nrows,
                         facet_scales = input$lesion_corr_facet_scales))
        }
      }



      if(length(args)>0){
        ## Use args to call the `plot_divset_cor` func
        divset_plot_obj <- do.call(plot_divset_cor, args)
        ##extract plot and dataframes from returned list
        p <- divset_plot_obj$plot
        df <- divset_plot_obj$exp_data
        divset_lesion_data <- divset_plot_obj$lesion_data
        if(troubleshoot_prints){ print(args)}
      }








      # actions to be taken if experiment == time series exp when button clicked
      }else if(input$experiment == "Time-Series Expression"){


        ##set argument default values to use if unless ortherwise specified by input
        input_plot_type <- 'line'
        single_panel = FALSE
        input_facet_nrow =2
        input_facet_scales = 'free_y'
        up_down_input = 'both'

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


        #modify DEGs input if specified
        if(isTruthy(input$degs_input)){
          if(grepl('upregulated',input$degs_input)){
            up_down_input = 'up'
            }else if(grepl('downregulated',input$degs_input)){
              up_down_input = 'down'
            }
          }

        #modify args for plot type
        if(isTruthy(input$time_series_plot_type)){
          if(input$time_series_plot_type =="Heatmap"){
            input_plot_type = 'heatmap'
            }else if(input$time_series_plot_type=="Single Panel Line Plot"){
              single_panel = TRUE
            }else if(input$time_series_plot_type=="Multi Panel Line Plot"){

              #Set multi-panel line plot specific args
              input_facet_nrow = as.integer(input$time_series_facet_nrows)
              input_facet_scales = ifelse(input$timeseries_facet_scales=='Gene-specific Y-axis',
                                      'free_y', 'free_x')
            }
        }
        if(troubleshoot_prints){ print('setting time-series args')}

        ## define args for time-series plot call
        args <- c(input_list(), #input_list of gene selection criteria
                list(
                #max number of genes to return
                max_n = as.integer(input$max_n),
                #return expr in response to which fungi
                fungi = input_fungi,
                #is mock expression required?
                include_mock = input$include_mock =='Yes',
                #DEGs only or all genes
                overlap_DEGs_only = input$degs_input !="No, show all genes",
                #If just DEGs, do we filter by direction of DE?
                up_down_only = up_down_input,
                #what type of plot? heatmap, line plot
                plot_type = input_plot_type,
                #label size
                strip.text.x_size = input$time_series_label_size,
                #if line plot, single panel or multipanel (facet_wrap)
                single_panel = single_panel,
                #if multipanel, how many rows?
                facet_nrow=input_facet_nrow,
                #if multipanel, should scale be consisent or gene-specific
                facet_scales=input_facet_scales)
                )


      if(length(args)>0){
        #call `plot_timeseries_expr` func
        timeseries_plot_obj <- do.call(plot_timeseries_expr, args)
        #extract data and plot from returned list
        df <- timeseries_plot_obj$data
        p <- timeseries_plot_obj$plot

        if(troubleshoot_prints){ print(args)}
      }



    ## actions after generate results button if experiment chosen is GRN
    }else if (input$experiment == "Gene Regulatory Network Analysis") {
      ##if blank gene input, return df of all hub stats
      if(trimws(input$geneInput)==""){
        grn_df <- db_query("SELECT TF, n_targets as TotalOutdegrees,
                        CASE
                       WHEN (At_ShortName IS NULL OR At_ShortName LIKE 'AT%G%') THEN TF
                       ELSE 'Ls' || At_ShortName
                       END AS name
                       FROM grn_hubs WHERE TotalOutdegrees>=40")
      }else{
        ## define arguments to get GRN table
        args <- c(input_list(), ## input for gene selection criteria
                  list(n_TFs = input$n_TFs, ## number of
                       min_subset_targets = input$min_targets,
                       order_by_column =  ifelse(input$grn_output_type=='Individual edges',
                                                 'SumImportance',input$hubs_sortby),
                       return_edges = ifelse(input$grn_output_type=='Individual edges',TRUE,FALSE)))
        grn_df <- do.call(get_pred_regulators, args)

        n_hubs <- length(unique(grn_df$TF))
      }


      #################################################################
      ### Make a legend for GRN table based on GRN output types   #####
      #################################################################


      if (input$grn_output_type == 'Individual edges') {
        text <- paste0(
          "Table containing predicted gene regulatory network (GRN) individual edges for the top ",n_hubs,
          "regulators (by SumImportance) of your selected gene set, where:<br>",
          "<ul>",
          "<li><b>TF</b> and <b>Target</b>: Represent lettuce geneIDs of the predicted GRN edge.</li>",
          "<li><b>TFName</b> and <b>TargetName</b>: Denote shortnames based on the predicted Arabidopsis orthologue of the lettuce gene for the transcription factor (TF) and target, respectively.</li>",
          "<li><b>Importance</b>: Indicates the predicted confidence of the transcriptional interaction.</li>",
          "<li><b>wigwam_module_TF</b> and <b>wigwam_module_Target</b>: Time-series co-expression modules that consider expression post fungal infection.</li>",
          "<li><b>Ss_DivSet_coexp_cor</b> and <b>Bc_DivSet_coexp_cor</b>: Capture the Pearson's correlation coefficient of co-expression between the TF and Target in a diversity panel after <i>Sclerotinina sclerotiorum</i> and <i>Botrytis cinerea</i> infection, respectively.</li>",
          "</ul>",
          "For comprehensive details on the GRN and co-expression modules, refer to <a href='https://doi.org/10.1101/2023.07.19.549542'>Pink et al (2023)</a>. For further information on diversity panel expression, see <a href='https://doi.org/10.1007/s00122-022-04129-5'>Pink et al (2022)</a>."
        )

      }else if (input$grn_output_type == 'Aggregated regulator statistics') {
        text <- paste0(
          "The table showcases the top ",n_hubs ,
          " regulators for the genes based on user-selected criteria, which regulate a minimum of ",
          input$min_targets, " within the chosen genes and are sorted by ", input$hubs_sortby,
          ". Specifically:<br>",
          "<ul>",
          "<li><b>TF</b>: Represents the lettuce geneID.</li>",
          "<li><b>TFName</b>: Provides the shortname derived from the closest Arabidopsis orthologue.</li>",
          "<li><b>TotalOutdegrees</b>: Denotes the count of genes anticipated to be governed by this TF throughout the entire GRN, without limiting to a gene subset.</li>",
          "<li><b>SubsetOutdegrees</b>: Indicates the quantity of genes in the user-specified gene subset.</li>",
          "<li><b>ExpectedSubsetOutdegrees</b>: Represents the predicted number of genes this TF might regulate within the subset based on random probability.</li>",
          "<li><b>SubsetProportion</b>: The fraction of genes from the user-picked subset regulated by the TF.</li>",
          "<li><b>TFProportion</b>: The fraction of genes influenced by the TF (across the entire GRN) that are part of this subset.</li>",
          "<li><b>SumImportance</b>: Aggregates the confidence scores (importance) of interactions for all edges directed by the TF in this gene subset.</li>",
          "<li><b>AvgImportance</b>: Computes the mean importance across these edges.</li>",
          "<li><b>FoldEnrichment</b>: Calculates how many times more edges within the gene subset are observed in comparison to random chance, specifically defined as SubsetOutdegrees divided by ExpectedSubsetOutdegrees.</li>",
          "</ul>",
          "For comprehensive details on the gene regulatory network, refer to <a href='https://doi.org/10.1101/2023.07.19.549542'>Pink et al (2023)</a>.")
      }
    }

    if(troubleshoot_prints){ print('assigning plot to reactive variable')}

    ## update reactive values

    make_plot(p)
    plot_data(df)
    lesion_size_data(divset_lesion_data)
    tableData(grn_df)
    LegendText(text)
    if(troubleshoot_prints){ print('reactive values assigned')}
  })




  ######################################################################
  ####  Render Outputs - Plots, Tables, Legends of Data downloads   ####
  ######################################################################

  #these outputs are set as conditional panels in the UI
  #so will only be displayed when user selects a given experiment
  #but would still throw errors tried to render a value that didn't exist (even if not displayed in UI)
  #thats why we initilised all reactive values as NULL




  ## render lesion corr plot

  output$lesionCorrPlot <- renderPlot({
    make_plot()
  }, height = 710, width = function() { return(session$clientData$output_lesionCorrPlot_width) })



  ## render timeseries plot

  output$timeSeriesPlot <- renderPlot({
    make_plot()
  },
  height = 710,
  width = function() { return(session$clientData$output_timeSeriesPlot_width) })


  # Render the table based on tableData reactiveVal
  output$mainTable <- renderDT({
    datatable(tableData(),
              extensions = 'Buttons',
              options = list(
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
              )
    )
  })


  ## render table legend
  output$dynamicLegend <- renderUI({
    return( HTML(LegendText()))
  })


  ## download handler to manage download of plot expression data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(gsub(' ','-',tolower(input$experiment)),'-data-', Sys.Date(), '.csv', sep='')
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


  ## Modal (pop-up) to manage options (dimesions, resolution, etc...) to download plots

  observeEvent(input$showDownloadOptions, {
    # Show modal dialog when the actionButton is clicked
    showModal(modalDialog(
      title = "Download Options",
      #plot file extension choice...
      radioButtons("fileType", "Choose File Type", choices = c("png", "tiff", "pdf", "RDATA")),
      numericInput("imgWidth", "Image Width", value = 1000, min = 300, max = 5000),
      numericInput("imgHeight", "Image Height", value = 750, min = 200, max = 5000),
      numericInput("imgResolution", "Resolution (DPI)", value = 100, min = 72, max = 1200),
      footer = tagList(
        downloadButton(outputId = "save_plot", "Save Plot"),
        modalButton("Cancel")
      ),
      easyClose = TRUE
    ))
  })


  ## Download
  output$save_plot <- downloadHandler(
    filename = function() {
      # Store user's modal choices in the reactiveValues
      plotDownloadUserChoices$fileType <- input$fileType
      plotDownloadUserChoices$imgWidth <- input$imgWidth
      plotDownloadUserChoices$imgHeight <- input$imgHeight
      plotDownloadUserChoices$imgResolution <- input$imgResolution

      ##construct default filename
      fname <- paste(gsub(' ','-',tolower(input$experiment)),"-plot-", Sys.Date(), ".", plotDownloadUserChoices$fileType, sep = "")
      if(troubleshoot_prints){
        print(paste("Filename:", fname))
        }
      return(fname)
    },
    ## content of plot download that handles several file extensions
    content = function(file) {
      if (plotDownloadUserChoices$fileType == "png") {
        png(file, width = plotDownloadUserChoices$imgWidth, height = plotDownloadUserChoices$imgHeight, res = plotDownloadUserChoices$imgResolution)
        print(make_plot())
        dev.off()
      } else if (plotDownloadUserChoices$fileType == "tiff") {
        tiff(file, width = plotDownloadUserChoices$imgWidth, height = plotDownloadUserChoices$imgHeight, res = plotDownloadUserChoices$imgResolution)
        print(make_plot())
        dev.off()
      } else if (plotDownloadUserChoices$fileType == "pdf") {
        pdf(file, width = plotDownloadUserChoices$imgWidth/100, height = plotDownloadUserChoices$imgHeight/100)  # PDF uses inches
        print(make_plot())
        dev.off()
      } else if (plotDownloadUserChoices$fileType == "RDATA") {
        saveRDS(make_plot(), file)
      }

      # Close the modal after the user chooses options
      removeModal()
    }
  )








}
