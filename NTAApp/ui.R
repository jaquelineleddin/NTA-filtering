library(shiny)
library(bslib)


loadPage <- page_sidebar(
 sidebar = sidebar(width = 300, open = "always",
      fileInput("file", "File", buttonLabel = "Load path", accept = c("text/csv", ".txt")),
      splitLayout(cellWidths = c("60%", "40%"),
        radioButtons("sep", "Separator", choices = c("tab", ";", ",")),
        numericInput("skipNum", "Rows to skip", 4, min = 0)
        ),
      actionButton("loadTable", "Load Data Frame"),
      numericInput("rows", "Number of rows in preview", 10, min = 1),
      checkboxInput("deleteNegCCS", "Remove features with -1 for mobility"),
      checkboxInput("deleteDefault", "Remove columns not needed for filtering"),
      selectInput("deleteCol", "Manual remove of columns", c(), multiple = TRUE),
    ),
 card(
   card_header("Current data frame"),
      tableOutput("table")
    )
  )

parameterPage <- page_fluid(
  card(
    card_header("Samples and identifyer"),
    splitLayout(
      selectInput("firstSampleCol", "First column of all samples (incl. QC and blanks)",
                  choices = c(), size = 10, selectize = FALSE),
      numericInput("numOfSamples", "Number of all samples (incl. QC and blanks)", 1, min = 1)
      ),
    splitLayout(
      textInput("ID.QC", "ID of quality control", placeholder = "QC"),
      textInput("ID.blanks", "ID of blanks", placeholder = "Blank")
      ),

      verbatimTextOutput("qc_columns"),
      verbatimTextOutput("blank_columns"),
      verbatimTextOutput("sample_columns")

    ),
  card(
    card_header("Missing values"),
    helpText("In data preprocessed by MSDail missing values in sample columns are marked as zeros.
             For further processing they must be replaced by NA. Same for negative values."),
    actionButton("repNegValuesBtn", "Replace", width = "20%"),
    tableOutput("replaced.table")
  )
)

plotOverviewPage <- navset_card_tab(
  sidebar = sidebar(width = 250, open = "always",
                    title = "Loading overview plots...please wait!",
                    numericInput("maxIntensityFilter", "Lower threshold for feature intensity", 10000, min = 0),
                    numericInput("maxNumberZeroFilter", "Maximum number of missing values in QC", 1, min = 0),
                    numericInput("minRSDFilter", "Maximum of RSD% in QC", 35, min = 1),
                    actionButton("filterBtn", "Filter rows"),
                    textOutput("filterMsg"),
                    downloadButton("savedata.afterFilter", "Save filtered data"),
                    downloadButton("save.overviewplot", "Save overview plot")
                    ),
  nav_panel(
    title = "Before",
    card(
      card_header("Maximum intensities"),
      plotOutput("maxIntPlot")
    ),
    card(
      card_header("Missing values in quality control"),
      plotOutput("qcZerosPlot")
    ),
    card(
      card_header("Relative standard deviation (RSD%) in quality control"),
      plotOutput("rsdPlot")
    )
  ),
  nav_panel(
    title = "After",
    card(
      card_header("Maximum intensities"),
      plotOutput("maxIntPlot.filtered")
    ),
    card(
      card_header("Missing values in quality control"),
      plotOutput("qcZerosPlot.filtered")
    ),
    card(
      card_header("Relative standard deviation (RSD%) in quality control"),
      plotOutput("rsdPlot.filtered")
    )
  )
)

boolFilterPage <- page_sidebar(
  sidebar = sidebar(width = 300, open = "always",
                    radioButtons("zerosSamplesChoice", "Count missing values in samples?",
                                 choices = c("No", "Yes")),
                    uiOutput("zeroSamplesUI"),
                    selectInput("matchFilterSlct",
                                "Select columns to filter out features, if the match is false",
                                choice = c(), multiple = TRUE),
                    numericInput("rowNum.filtered", "Number of rows in preview", 10, min=1),
                    textOutput("filterMsg.optfilter"),
                    downloadButton("save.afterBoolFilter", "Save current data frame")
                    ),
  uiOutput("zeroSamplesPlotUI"),
  card(
    card_header("Filtered data frame"),
    tableOutput("data.filtered")
  )
)

normalizePage <- page_sidebar(
  sidebar = sidebar(width = 300, open = "always",
                    radioButtons("normalizeChoice", "Normalize intensities of samples?",
                                 choices = c("Yes", "No"), selected = "No"),
                    uiOutput("normalizeUI")),
  uiOutput("normalizeCards")
)

namesPage <- page_sidebar(
  sidebar = sidebar(width = 300, open = "always",
                    selectInput("newNamesSlct", "Columns to combine for new feature names",
                                multiple = TRUE, choices = c()),
                    numericInput("rowNum.newnames", "Number of rows in preview", 10, min=1),
                    downloadButton("save.newnamesdf", "Save data frame")
                    ),
  card(
    card_header("Current data frame"),
    tableOutput("newnamesdf")
  )
)

transposePage <- page_sidebar(
  sidebar = sidebar(width = 300, open = "always",
                    selectInput("transpose.rowname", "Column to use as name of features",
                                choices = c(),  size = 10, selectize = FALSE),
                    textInput("classes", "Class names", placeholder = "QC, A, B, C"),
                    radioButtons("qcInClasses", "With QC?", choices = c("Yes", "No"), selected = "No"),
                    radioButtons("blanksInClasses", "With Blanks?", choices = c("Yes", "No"), selected = "No"),
                    actionButton("transposeBtn", "Transpose"),
                    downloadButton("save.transposeddf", "Save data frame")
                    ),
  card(
    card_header("Transposed data frame"),
    tableOutput("transposeddf")
  )
)


page_navbar(
  id = "main",
  title = "NTA Workflow",
  nav_panel("Load Data",
           loadPage),
  nav_panel("Samples",
           parameterPage),
  nav_panel("Filter",
            plotOverviewPage),
  nav_panel("Optional filter",
            boolFilterPage),
  nav_panel("Normalization",
            normalizePage),
  nav_panel("Feature names",
            namesPage),
  nav_panel("Transpose",
            transposePage)
)

