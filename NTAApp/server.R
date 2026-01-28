library(shiny)
library(NTAWorkflow)
library(ggplot2)
library(gridExtra)

function(session, input, output){

  nav_hide("main", "Samples")
  nav_hide("main", "Filter")
  nav_hide("main", "Optional filter")
  nav_hide("main", "Normalization")
  nav_hide("main", "Feature names")
  nav_hide("main", "Transpose")

  # change here maximum file size of source data frame (default is 70)
  options(shiny.maxRequestSize = 70*1024^2)

  # load table tab
  df <- eventReactive(input$loadTable,
                      {
                        req(input$file)

                        if(input$sep == "tab") sep = "\t" else sep = input$sep

                        data <- NTAWorkflow::load_data(input$file$datapath,
                                                       skip = input$skipNum,
                                                       sep = sep,
                                                       drop_last_cols = 0,
                                                       drop_cols = NULL)

                        updateSelectInput(session, "deleteCol", choices = colnames(data))
                        updateNumericInput(session, "rows", max = nrow(data))

                        return(data)
                      })

  current.df <- reactive({
    req(df())

    colsToRemove = c("Post.curation.result", "Fill..", "Reference.RT",
                     "Reference.m.z.", "Formula", "Ontology", "INCHIKEY",
                     "SMILES", "Annotation.tag..VS1.0.", "Comment",
                     "Manually.modified.for.quantification",
                     "Manually.modified.for.annotation",	"Isotope.tracking.parent.ID",
                     "Isotope.tracking.weight.number",	"Total.score",
                     "RT.similarity", "Dot.product",	"Reverse.dot.product",
                     "Fragment.presence..", "S.N.average",
                     "Spectrum.reference.file.name",	"MS1.isotopic.spectrum",
                     "MS.MS.spectrum", "X1", "X2")
    
    if(input$deleteDefault){
      allColsToRemove <- union(input$deleteCol, colsToRemove)
      updateSelectInput(session, "deleteCol", selected = allColsToRemove)
    }

    if(is.null(input$deleteCol)){
      data <- df()
      cols <- colnames(df())}
    else{
      col.ind = !(colnames(df()) %in% input$deleteCol)
      data <- df()[, col.ind]
      cols <- colnames(df())[col.ind]
    }
    
    if(input$deleteNegCCS){
      data <- remove_neg_mobility(data)
    }
    
    updateSelectInput(session, "firstSampleCol", choices = cols)

    return(data)
  })


  output$table <- renderTable(current.df()[1:input$rows,], striped = TRUE)


  # parameter - list of sample columns

  observeEvent(current.df(),{
    if(is.null(current.df()))
      nav_hide("main", "Samples")
    else
      nav_show("main", "Samples")

               })

  stringList <- function(str.list){
    res <- ""

    for(item in str.list){
      res <- paste(res, item, sep = "\n")
    }
    return(res)
  }

  all_samples <- reactive({
    req(input$firstSampleCol)

    first.sample.ind = which(colnames(current.df()) == input$firstSampleCol)
    colnames(current.df())[first.sample.ind:(first.sample.ind+input$numOfSamples-1)]
  })

  qc_samples <- reactive({
    req(all_samples())

    if(input$ID.QC == "")
      c("")
    else
      grep(input$ID.QC, all_samples(), value = TRUE)
  })

  blank_samples <- reactive({
    req(all_samples())

    if(input$ID.blanks == "")
      c("")
    else
      grep(input$ID.blanks, all_samples(), value = TRUE)
  })

  only_samples <- reactive({
    req(all_samples())

    res <- all_samples()
    res <- res[!(res %in% blank_samples())]
    res <- res[!(res %in% qc_samples())]

    return(res)
  })

  output$qc_columns <- renderText(paste0("Columns of quality control:\n", stringList(qc_samples())))
  output$blank_columns <- renderText(paste0("Columns of blank:\n", stringList(blank_samples())))
  output$sample_columns <- renderText(paste0("Columns of samples:\n", stringList(only_samples())))

  # Replace negative intensities and missing values
  replaced.df <- eventReactive(input$repNegValuesBtn,
                               {
                                 df <- NTAWorkflow::replace_missing_values(current.df(), all_samples())
                                 collogi <- c()
                                 colnewnames <- c()


                                 for(col in colnames(df)){
                                   if(class(df[[col]])=="logical"){
                                      collogi <- append(collogi, col)
                                   } else {
                                     colnewnames <- append(colnewnames, col)
                                   }
                                 }

                                 # remove all sample columns for list of new names
                                 colnewnames <- colnewnames[!(colnewnames %in% all_samples())]

                                  updateSelectInput(session, "matchFilterSlct", choices = collogi)
                                  updateSelectInput(session, "newNamesSlct", choices = colnewnames)

                                  return(df)
                               })


  output$replaced.table <- renderTable(replaced.df()[1:input$rows,], striped = TRUE)

  observeEvent(replaced.df(),{
    if(is.null(replaced.df()))
      nav_hide("main", "Filter")
    else
      nav_show("main", "Filter")

  })


  # Plots for data overview
  maxIntensity <- reactive({

    req(only_samples())
    req(replaced.df())

    NTAWorkflow::calc_max_intensity(replaced.df(), only_samples())
    })

  output$maxIntPlot <- renderPlot({
    gridExtra::grid.arrange(maxIntensity()$plot, maxIntensity()$subplot, ncol = 2)
    })

  qcZeros <- reactive({
    req(qc_samples())
    req(replaced.df())

    NTAWorkflow::count_missing_values(replaced.df(), qc_samples())
  })

  output$qcZerosPlot <- renderPlot(qcZeros()$plot)

  rsd <- reactive({
    req(qc_samples())
    req(replaced.df())

    NTAWorkflow::calc_RSD(replaced.df(), qc_samples())
  })

  output$rsdPlot <- renderPlot({
    gridExtra::grid.arrange(rsd()$plot, rsd()$subplot, ncol = 2)
  })

  beforeFilter.df <- reactive({
    req(replaced.df())

    cbind(replaced.df(),
          data.frame(
            "max.Intensity" = maxIntensity()$max_intensity,
            "missing.values" = qcZeros()$missing_values,
            "rsd" = rsd()$RSD))

  })

  # Filter data

  afterFilter.df <- eventReactive(input$filterBtn,
                                  {
                                    NTAWorkflow::filter_dataframe(beforeFilter.df(),
                                                                  colnames = c("max.Intensity",
                                                                               "missing.values",
                                                                               "rsd"),
                                                                  thresholds = c(input$maxIntensityFilter,
                                                                                 input$maxNumberZeroFilter,
                                                                                 input$minRSDFilter))
                                    })


  output$filterMsg <- renderText({
    req(beforeFilter.df())
    req(afterFilter.df())

    numOfdelRows <- nrow(beforeFilter.df()) - nrow(afterFilter.df())
    percentage <- round(numOfdelRows / nrow(beforeFilter.df()) * 100, digits = 0)

    paste("Deleted", numOfdelRows, "of", nrow(beforeFilter.df()),"(", percentage,"%)", "features.")
  })


  observeEvent(afterFilter.df(),{
    if(is.null(afterFilter.df())){
      nav_hide("main", "Optional filter")
      nav_hide("main", "Normalization")
      nav_hide("main", "Feature names")
      nav_hide("main", "Transpose")
    }
    else {
      nav_show("main", "Optional filter")
      nav_show("main", "Normalization")
      nav_show("main", "Feature names")
      nav_show("main", "Transpose")
    }
  })


  maxIntensity.filtered <- reactive({

    req(only_samples())
    req(afterFilter.df())

    NTAWorkflow::calc_max_intensity(afterFilter.df(), only_samples())
  })

  output$maxIntPlot.filtered <- renderPlot({
    maxIntensity.filtered()$plot
  })

  qcZeros.filtered <- reactive({
    req(qc_samples())
    req(afterFilter.df())

    NTAWorkflow::count_missing_values(afterFilter.df(), qc_samples())
  })

  output$qcZerosPlot.filtered <- renderPlot(qcZeros.filtered()$plot)

  rsd.filtered <- reactive({
    req(qc_samples())
    req(afterFilter.df())

    NTAWorkflow::calc_RSD(afterFilter.df(), qc_samples())
  })

  output$rsdPlot.filtered <- renderPlot({
    rsd.filtered()$plot
  })

  output$savedata.afterFilter <- downloadHandler(

    filename = function(){"data_afterfirstfilter.csv"},
    content = function(file){
      write.table(afterFilter.df(), file, sep="\t", row.names = FALSE, quote = FALSE, na="")
    }
  )

  output$save.overviewplot <- downloadHandler(

    filename = function(){"overview_beforefilter.png"},
    content = function(file){
      overview_plot = grid.arrange(maxIntensity()$plot, maxIntensity()$subplot,
                                   qcZeros()$plot,
                                   rsd()$plot + guides(fill = "none"),
                                   rsd()$subplot,
                                   layout_matrix = matrix(c(1,2,3,3,4,5),
                                                          nrow = 3,
                                                          ncol = 2,
                                                          byrow = TRUE))

      ggsave(filename = file, plot = overview_plot,
             device = "png", width = 12, height = 8)
    }
  )

  # create plot for filter by missing values in samples
  sampleZero <- reactive(
    {
      req(afterFilter.df())

      if(input$zerosSamplesChoice == "Yes")
        NTAWorkflow::count_missing_values(afterFilter.df(), only_samples())

    })

  output$zeroSamplesPlot <- renderPlot(sampleZero()$plot)

  # Optional filter by bool
  matchfiltered.df <- reactive({
    req(afterFilter.df())

    df <- afterFilter.df()

    if(input$zerosSamplesChoice == "Yes"){
      df$miss.values.samples <- sampleZero()$missing_values
      df <- NTAWorkflow::filter_dataframe(df,
                                        colnames = c("miss.values.samples"),
                                        thresholds = c(input$maxZerosSamples),
                                        modes = c("<="))
    }

    if(!is.null(input$matchFilterSlct))
      df <- NTAWorkflow::filter_matched(df,input$matchFilterSlct)

    return(df)
    })


  output$data.filtered <- renderTable(matchfiltered.df()[1:input$rowNum.filtered,], striped = TRUE)

  output$filterMsg.optfilter <- renderText({
    req(afterFilter.df())
    req(matchfiltered.df())

    numOfdelRows <- nrow(afterFilter.df()) - nrow(matchfiltered.df())
    percentage <- round(numOfdelRows / nrow(afterFilter.df()) * 100, digits = 0)

    paste("Deleted", numOfdelRows, "of", nrow(afterFilter.df()),"(", percentage,"%)", "features.")
  })

  # filter missing values in samples
  output$zeroSamplesUI <- renderUI({
    if(input$zerosSamplesChoice == "Yes"){
      numericInput("maxZerosSamples", "Maximum number of missing values in samples",
                   1, min = 0)
    }
  })

  output$zeroSamplesPlotUI <- renderUI({
    if(input$zerosSamplesChoice == "Yes"){
      card(
        card_header("Missing values in samples"),
        max_height = 300,
        plotOutput("zeroSamplesPlot")
      )
    }
  })

  # save df after bool filter
  output$save.afterBoolFilter <- downloadHandler(

    filename = function(){"data_aftermatchfilter.csv"},
    content = function(file){
      write.table(matchfiltered.df(), file, sep="\t", row.names = FALSE, quote = FALSE, na="")
    }
  )

  # optional: normalization of sample intensities
  output$normalizeUI <- renderUI({
    if(input$normalizeChoice == "Yes"){
        tagList(fileInput("norm.factor.file", "Load data frame (first column: sample name, second column: factor)",
                buttonLabel = "Load path", accept = c("text/csv", ".txt")),
                radioButtons("sep.factor", "Separator", choices = c("tab", ";", ",")),
                actionButton("load.norm.df.Btn", "Load data frame"),
                actionButton("normalize.Btn", "Normalize"),
                downloadButton("save.normdf", "Save current data frame"))
    }
  })

  output$normalizeCards <- renderUI({
    if(input$normalizeChoice == "Yes"){
      tagList(card(
        card_header("Factors for normalization"),
        tableOutput("norm.factors.table")
      ),
      card(
        card_header("Normalized data frame"),
        tableOutput("normalized.table")
      ))
    }
  })

  # load normalization data frame
  norm.factors <- eventReactive(input$load.norm.df.Btn,
                      {
                        req(input$norm.factor.file)

                        if(input$sep.factor == "tab") sep = "\t" else sep = input$sep.factor

                        read.csv(input$norm.factor.file$datapath, sep = sep)
                      })

  output$norm.factors.table <- renderTable(norm.factors(), striped = TRUE)

  normalized.df <- eventReactive(input$normalize.Btn,
                                 { req(norm.factors())
                                   req(matchfiltered.df())

                                   NTAWorkflow::normalize_intensity(matchfiltered.df(),
                                                                    norm.factors())
                                 })

  output$normalized.table <- renderTable(normalized.df(), striped = TRUE)

  # save data frame after normalization
  output$save.normdf <- downloadHandler(

    filename = function(){"data_normalized.csv"},
    content = function(file){
      write.table(normalized.df(), file, sep="\t", row.names = FALSE, quote = FALSE, na="")
    }
  )


  # add new names for features
  newnames.df <- reactive({
    req(matchfiltered.df())

    if(input$normalizeChoice == "Yes"){
      df <- normalized.df()
    } else
    {
      df <- matchfiltered.df()
    }

    # new names?
    if(length(input$newNamesSlct) > 1){
      # get list of column numbers
      nameslist <- c()

      for(col in input$newNamesSlct){
        colind = which(colnames(df) == col)
        nameslist <- append(nameslist, colind)
      }

      df <- NTAWorkflow::add_feature_names(df, nameslist)
    }

    # update names list for transpose
    currentCols = colnames(df)[!(colnames(df) %in% all_samples())]
    strCols = c()

    for(col in currentCols){
      if(class(df[[col]])=="character"){
        strCols <- append(strCols, col)
      }
    }

    updateSelectInput(session, "transpose.rowname", choices = strCols)

    return(df)
  })

  output$newnamesdf <- renderTable(newnames.df()[1:input$rowNum.newnames, ])

  output$save.newnamesdf <- downloadHandler(

    filename = function(){"data_newnames.csv"},
    content = function(file){
      write.table(newnames.df(), file, sep="\t", row.names = FALSE, quote = FALSE, na="")
    }
  )

  #transpose
  transpose.df <- eventReactive(input$transposeBtn,{

    req(newnames.df())

    # parse classes, remove white space at the beginning
    classes = gsub("^\\s+|\\s+$","",unlist(strsplit(input$classes, ",", fixed = TRUE)))

    samples <-  all_samples()

    if(input$qcInClasses == "No"){
      samples <- samples[!(samples %in% qc_samples())]
    }

    if(input$blanksInClasses == "No"){
      samples <-  samples[!(samples %in% blank_samples())]
    }

    df <- NTAWorkflow::transpose_df(newnames.df(),
                              colnum_of_name = which(colnames(newnames.df()) == input$transpose.rowname),
                              all_sample_cols = samples,
                              classes = classes)

    return(df)
  })

  output$transposeddf <- renderTable(transpose.df())


  # save data frame after transpose

  output$save.transposeddf <- downloadHandler(

    filename = function(){"data_transposed.csv"},
    content = function(file){
      write.table(transpose.df(), file, sep=",", row.names = FALSE, quote = TRUE, na="")
    }
  )






session$onSessionEnded(function() {
  stopApp()
})

}
