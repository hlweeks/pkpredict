#' Shiny application for pharmacokinetic prediction model
#'
#' Leave as-is? Or should any values be considered parameters
#' Currently specific to piperacillin (including some associated text)
#'
#' @return
#' Launches an interaction interface for implementing Bayesian pharmacokinetic modeling
#'
#' @export
#'
#' @examples

shiny_pkm <- function(){
  library(pkpredict)

  # Change the shiny app to only rely on functions/syntax within pkpredict package

  ui <- dashboardPage(
    dashboardHeader(
      titleWidth = 250,
      title = "Piperacillin TDM"
    ),
    dashboardSidebar(
      width = 250,
      tags$head(tags$style(HTML(".wrapper {overflow: visible !important;}"),
                           HTML(".sidebar { width: 240px; font-size: 18px;  margin-left: 5px; }"),
                           HTML(".main-header .logo { font-size: 30px; }"))),
      hr(),
      p("This application is designed to provide individualized estimates of drug exposure for
        critically ill patients with sepsis that have received one or more doses of piperacillin."),
      p("Enter patient dosing schedule (historical and/or proposed) and
        measurements of drug concentration in blood (if any)."),
      hr(),

      checkboxInput("common", "Common Dosing Pattern"),
      em("Manual entry of infusion data is disabled while this box is checked."),

      conditionalPanel(
        condition = "input.common == true",
        sliderInput("num", "Number of Doses", min = 1, max = 10, value = 5, step = 1),
        numericInput("freq", "Dose Frequency (h)", value = 8, min = 0),
        numericInput("duration", "Infusion Duration (h)", value = .5, min = 0),
        numericInput("infRate", "Infusion Rate (g/h)", value = 6, min = 0)
      )

      ),
    dashboardBody(
      tags$head(tags$style(HTML("#thres { font-size: 20px; height: 150%; }"),
                           HTML("#dosing { font-size: 20px; }"),
                           HTML("#sample { font-size: 20px; }"))),
      tags$style(type = "text/css", "label { font-size: 18px; }"),

      h3("Piperacillin Therapeutic Drug Monitoring"),

      fluidRow(
        column(width = 6,
               h4("Infusion Schedule", style = "font-size: 20px;"),
               rHandsontableOutput("dosing")
        ),
        column(width = 6,
               h4("Concentration Data", style = "font-size: 20px;"),
               rHandsontableOutput("sample")
        )
      ),

      HTML('<br/>'),
      checkboxInput("MCMC", "Compute credible interval for fraction of time above threshold statistic using Markov chain Monte Carlo sampling of the posterior distribution.",
                    width = "200%"),
      actionButton("goPlot", "Update Plot", style = 'padding: 10px; font-size: 20px'),
      HTML('<br/><br/>'),

      numericInput('thres', "Threshold (μg/ml)", value = 64, min = 0, width = '175px'),
      em("Threshold value is typically some multiple of the minimum inhibitory concentration
         (MIC) for the target microorganism.", style = "font-size: 18px;"),
      HTML('<br/><br/>'),
      em("'fT > threshold' in the plot legend provides an estimate of the fraction of time
         spent above the specified threshold.", style = "font-size: 18px;"),

      plotOutput("plot", hover = "plot_hover"),
      verbatimTextOutput("info"),

      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }")
      )
    )


  lpk_mean_d <- c(lv_1=3.223, lk_10=-1.650, lk_12 = -7, lk_21 = -7); ler_mean_d <- 2.33
  lpr_mean_d <- c(lpk_mean_d, ler_mean_d)

  # Define server logic required to draw plot
  server <- function(input, output) {
    # App title with line break
    output$titleText <- renderUI(HTML("Piperacillin Therapeutic Drug Monitoring"))


    #rhandsontable for dosing information
    output$dosing <- renderRHandsontable({
      if(1 - input$common){
        if(is.null(input$dosing)){
          doseDF  <- data.frame(
            "Start (h)" = c(0, 8, 16, 24, 32, rep(0,5)),
            "Duration (h)" = c(rep(0.5, 5), rep(0,5)),
            "Rate (g/h)" = c(6, 6, 6, 6, 6, rep(0,5)), check.names=FALSE)
        }else{
          doseDF <- hot_to_r(input$dosing)
        }
        rhandsontable(doseDF, colWidths = c(85,115,100), rowHeights = rep(30, 10),
                      colHeaders = c("Start (h)", "Duration (h)", "Rate (g/h)"))
      }else{
        comPat <- hot_to_r(input$dosing)
        nDose <- as.numeric(input$num)
        begin <- seq(0, (nDose - 1) * input$freq, by = input$freq)
        dur <- rep(input$duration, input$num)
        kR <- rep(input$infRate, nDose)
        comPat[1:input$num, "Start (h)"] <- begin
        comPat[1:input$num, "Duration (h)"] <- dur
        comPat[1:input$num, "Rate (g/h)"] <- kR
        if(input$num < 10){
          comPat[(input$num + 1):10,] <- matrix(0, ncol = 3, nrow = 10 - input$num)
        }

        rhandsontable(comPat, colWidths = c(85,115,100), rowHeights = rep(30, 10),
                      colHeaders = c("Start (h)", "Duration (h)", "Rate (g/h)"))
      }
    })

    #rhandsontable for sample information
    output$sample <- renderRHandsontable({
      vec <- numeric(10)
      sampDF = data.frame("Time (h)" = vec,
                          "Conc. (μg/ml)" = vec,
                          check.names=FALSE)
      rhandsontable(sampDF, colWidths = c(95, 135), rowHeights = rep(30, 10),
                    colHeaders = c("Time (h)", "Conc. (μg/ml)"))
    })


    ##########################################################

    # These functions allow me to grab the data from rhandsontable to use in plotting
    sampTable <- eventReactive(input$goPlot, {
      live_data = hot_to_r(input$sample)
      return(live_data)
    })

    doseTable <- eventReactive(input$goPlot, {
      live_data = hot_to_r(input$dosing)
      return(live_data)
    })


    ##########################################################

    # Create the plot
    output$plot <- renderPlot({
      if(input$goPlot == 0){
        dat <- data.frame("time_h" = numeric(0), "conc_mg_dl" = numeric(0))

        # DOSING INFORMATION
        ivtData <- list(list(begin=0.0, end=0.5, k_R=6),
                        list(begin=8.0, end=8.5, k_R=6),
                        list(begin=16.0, end=16.5, k_R=6),
                        list(begin=24.0, end=24.5, k_R=6),
                        list(begin=32.0, end=32.5, k_R=6))


        plot(pkm(conc_mg_dl ~ time_h, data = dat, ivt = ivtData))

      }else{

        # Get data from rHandsonTable
        stab <- sampTable()
        dtab <- doseTable()

        # Won't produce plot unless both sample information and dosing schedule has been provided
        if(sum(stab > 0) > 0 & sum(dtab > 0) > 0){

          # SAMPLE INFORMATION
          datHot <- stab[apply(stab, MARGIN = 1, function(x) any(x > 0)),]
          # Required for compatibility with functions
          names(datHot) <- c("time_h", "conc_mg_dl")

          # DOSING INFORMATION
          # Get data from rhandsontable
          ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
          # Convert duration of infusions to start/end times
          ivtHot[,"Duration (h)"] <- ivtHot[, "Start (h)"] + ivtHot[, "Duration (h)"]

          # Required for compatibility with functions
          names(ivtHot) <- c("begin", "end", "k_R")
          ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))
          ivtData <- ivtHot


          pkm_mod <- pkm(conc_mg_dl ~ time_h, data = datHot, ivt = ivtHot,
                         mcmc = isolate(input$MCMC), shiny = TRUE)

          plot(pkm_mod)

        }else if(sum(stab > 0) == 0){
          if(sum(dtab > 0) > 0){
            # DOSING INFORMATION
            # Get data from rhandsontable
            ivtHot <- dtab[apply(dtab, MARGIN = 1, function(x) any(x > 0)),]
            # Convert duration of infusions to start/end times
            ivtHot[,"Duration (h)"] <- ivtHot[, "Start (h)"] + ivtHot[, "Duration (h)"]

            # Required for compatibility with functions from Bayes.R
            names(ivtHot) <- c("begin", "end", "k_R")
            ivtHot <- apply(ivtHot, MARGIN = 1, function(x) list(begin = x[1], end = x[2], k_R = x[3]))
            ivtData <- ivtHot

            pkm_mod <- pkm(ivt = ivtHot, mcmc = isolate(input$MCMC), shiny = TRUE)
            plot(pkm_mod)

          }
        }
      }
      #Display coordinates when hovering over a point
      output$info <- renderText({
        paste("Time=", round(input$plot_hover$x,2), "h",
              "\nConcentration=", round(input$plot_hover$y,2), "μg/ml", sep=" ")
      })
    })

  }
  shiny::shinyApp(ui = ui, server = server)
}



