library(shiny)
library(DT)

funcs_path <-
source()


ui <- fluidPage(
  titlePanel("Integrated Analysis App"),

  sidebarLayout(
    sidebarPanel(
      selectInput("experiment", "Choose an Experiment",
                  choices = c("Lesion Size Correlation", "Time-Series Expression", "Gene Regulatory Network Analysis")),

      selectInput("geneSelectionMethod", "Gene selection method",
                  choices = c("Lettuce GeneID", "Ortholog of Arabidopsis Gene", "Genes with GO-term", "Genes with Protein Domain")),

      # Dynamic Default for Text Input based on Gene Selection Method
      uiOutput("geneInputUI"),

      # Dynamic UI elements based on experiment choice
      uiOutput("experimentFilters")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("mainPlot")),
        tabPanel("Table", DTOutput("mainTable"))
      )
    )
  )
)

server <- function(input, output, session) {

  # Render the appropriate default text based on gene selection method
  output$geneInputUI <- renderUI({
    switch(input$geneSelectionMethod,
           "Lettuce GeneID" = textInput("geneInput", "Enter Gene ID", value="LettuceGene12345"),
           "Ortholog of Arabidopsis Gene" = textInput("geneInput", "Enter Ortholog Gene ID", value="ArabidopsisGene12345"),
           "Genes with GO-term" = textInput("geneInput", "Enter GO Term", value="GO:0000001"),
           "Genes with Protein Domain" = textInput("geneInput", "Enter Protein Domain ID", value="PD:12345")
    )
  })

  output$experimentFilters <- renderUI({
    if (input$experiment == "Lesion Size Correlation") {
      selectInput("geneFilter1", "Filter by Gene", choices = c("GeneA", "GeneB"))

    } else if (input$experiment == "Time-Series Expression") {
      selectInput("geneFilter2", "Filter by Time", choices = c("Time1", "Time2"))

    } else if (input$experiment == "Gene Regulatory Network Analysis") {
      selectInput("geneFilter3", "Filter by Regulator", choices = c("Regulator1", "Regulator2"))
    }
  })

  output$mainPlot <- renderPlot({
    if (input$experiment == "Lesion Size Correlation") {
      plot(1:10, rnorm(10), main="Scatter Plot for Lesion Size Correlation")
    } else if (input$experiment == "Time-Series Expression") {
      plot(sin(seq(1, 10, by = 0.1)), type = "l", main="Time-Series Expression")
    } else {
      NULL
    }
  })

  output$mainTable <- renderDT({
    if (input$experiment == "Gene Regulatory Network Analysis") {
      datatable(data.frame(Gene=c("GeneA", "GeneB"), Regulator=c("Regulator1", "Regulator2")))
    } else {
      return(NULL)
    }
  })
}

shinyApp(ui = ui, server = server)
