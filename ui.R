# ui.R
library(shiny)
library(ggplot2)
library(plyr)
library(FD)

pageWithSidebar(
  headerPanel("traitR"),
  sidebarPanel(
    fileInput(
      'file1',
      'Choose CSV File',
      accept = c('text/csv', 'text/comma-separated-values,text/plain', '.csv')
    ),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE),
    
    radioButtons('sep', 'Separator',
                 c(
                   Comma = ',',
                   Semicolon = ';',
                   Tab = '\t'
                 ), ','),
    radioButtons(
      'quote',
      'Quote',
      c(
        None = '',
        'Double Quote' = '"',
        'Single Quote' = "'"
      ),
      '"'
    ),
    fluidRow(
      column(6, selectInput(
        "speciesCol", "Species Column", c("1" = "1")
      )),
      column(6, selectInput(
        "ancestralGrp", "Ancestral Species", c("1" = "1")
      )),
      column(6, selectInput(
        "attackedCol", "Attacked Column", c("1" = "1", "2" = "2")
      )),
      column(6, selectInput(
        "attackedGrp", "Attacked Identifier", c("1" = "1")
      )),
      column(6, selectInput(
        "subsetCol", "Column for Subsetting", c("1" = "1")
      )),
      column(6, selectInput(
        "subsetGrp", "Subset Identifier", c("1" = "1", "2" = "2")
      )),
      column(
        6,
        selectInput(
          "traitGrp",
          "Traits:",
          choices = NULL,
          selected = NULL,
          multiple = T,
          selectize = T
        )
      ),
      column(
        6,
        selectInput(
          "traitFac",
          "Factor Traits:",
          choices = NULL,
          selected = NULL,
          multiple = T,
          selectize = TRUE
        )
      ),
      column(
        6,
        selectInput(
          "traitOrdFac",
          "Ordered Factor Traits:",
          choices = NULL,
          selected = NULL,
          multiple = T,
          selectize = TRUE
        )
      ),
      column(
        6,
        numericInput(
          "niter",
          label = h4("Number of permutations:"),
          value = 1000
        ),
        h5("Perform traitR analysis"),
        wellPanel(actionButton("do", "Submit")),
        h5("Download traitR results"),
        wellPanel(downloadButton("downloadData", "Download"))
      )
    )
  ),
  mainPanel(tabsetPanel(
    tabPanel(
      "Plots and Analysis",
      plotOutput("plot1"),
      br(),
      tableOutput("sum"),
      br(),
      tableOutput("pvalueLoc"),
      br(),
      tableOutput("pvalueSpread"),
      h5(textOutput("speciescentroids")),
      tableOutput("centroidtable"),
      h5(textOutput("attackedcentroids")),
      tableOutput("centroidtable2"),
      plotOutput("plot2")
    ),
    tabPanel("PCA Output",
             tableOutput("filetable"))
  ))
)